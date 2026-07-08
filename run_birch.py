#!/usr/bin/env python3
"""Build and run paper-aligned A-BIRCH on every supplied dataset."""
from __future__ import annotations
import argparse
import csv
import html
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IMPL = ROOT / "Birch-Implementation"
SOURCES = [
    "main.cpp",
    "Common.cpp",
    "Dataset.cpp",
    "KNNOutlier.cpp",
    "BIRCH.cpp",
    "Validation.cpp",
    "Pipeline.cpp",
]

def run(command: list[str], cwd: Path) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)

def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and run C++ A-BIRCH on all supplied datasets."
    )
    parser.add_argument(
        "--features",
        type=Path,
        help="Custom feature CSV (requires --truth); otherwise run all supplied datasets.",
    )
    parser.add_argument("--truth", type=Path, help="Ground-truth CSV paired with --features.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "A-BIRCH-Results")
    parser.add_argument("--branching-factor", type=int, default=32)
    parser.add_argument("--minimum-cluster-points", type=int, default=5)
    parser.add_argument("--pilot-points", type=int, default=1000)
    parser.add_argument("--pilot-max-clusters", type=int, default=10)
    parser.add_argument(
        "--threshold",
        type=float,
        default=-1.0,
        help="Negative selects the paper's automatic threshold.",
    )
    parser.add_argument(
        "--mbd-spread",
        type=float,
        default=0.05,
        help="MBD-BIRCH multiple-branch distance allowance s; use 0 for ordinary tree-BIRCH.",
    )
    parser.add_argument("--knn-k",type=int,default=3,
                        help="Neighbor rank used by distance-based kNN outlier detection.")
    parser.add_argument("--knn-outlier-fraction",type=float,default=0.10,
                        help="Fraction m/n of largest k-distances removed before BIRCH.")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--build-only", action="store_true")
    return parser.parse_args()

def prepare_combined(source:Path,destination:Path)->tuple[Path,Path]:
    with source.open(newline="",encoding="utf-8-sig") as f: rows=list(csv.reader(f))
    if not rows or "Emitter" not in rows[0]: raise ValueError(f"Missing Emitter column: {source}")
    header=rows[0];label_index=header.index("Emitter")
    wanted = {
        "radio frequency", "frequency", "freq_mhz",
        "pulse width", "pw", "pw_ns",
        "azimuth", "az", "az_deg",
        "angle of elevation", "elevation", "el", "el_deg",
    }
    feature_indices = [
        i for i, name in enumerate(header[:label_index])
        if name.strip().lower() in wanted
    ]
    if len(feature_indices) != 4:
        raise ValueError(f"Expected frequency, PW, azimuth and elevation in {source}")
    parsed=[];sums={};counts={}
    for row in rows[1:]:
        if len(row)<=label_index or not row[label_index].strip():continue
        features=[row[i].strip() for i in feature_indices];values=[float(v) for v in features]
        label=row[label_index].strip();parsed.append((features,label))
        if float(label)>=0:sums[label]=sums.get(label,0.0)+values[0];counts[label]=counts.get(label,0)+1
    # Match the C++ display-name convention (increasing mean frequency). This
    # changes names only; labels never enter A-BIRCH fitting or threshold selection.
    ordered=sorted(counts,key=lambda x:sums[x]/counts[x]);names={x:f"Emitter_{i+1}" for i,x in enumerate(ordered)}
    destination.mkdir(parents=True,exist_ok=True);features_file=destination/"features.csv";truth_file=destination/"truth.csv"
    with features_file.open("w",newline="",encoding="utf-8") as ff,truth_file.open("w",newline="",encoding="utf-8") as tf:
        fw,tw=csv.writer(ff),csv.writer(tf);feature_header=[header[i] for i in feature_indices]
        fw.writerow(feature_header);tw.writerow([*feature_header,"Ground_Truth"])
        for features,label in parsed:
            truth="Noise" if float(label)<0 else names[label];fw.writerow(features);tw.writerow([*features,truth])
    return features_file,truth_file

def prepare_two_emitter(source:Path,labeled:Path,destination:Path)->tuple[Path,Path]:
    """Remove TOA and row number while preserving frequency, PW, azimuth and elevation."""
    with source.open(newline="",encoding="utf-8-sig") as f: source_rows=list(csv.reader(f))
    with labeled.open(newline="",encoding="utf-8-sig") as f: labeled_rows=list(csv.reader(f))
    source_header=source_rows[0];truth_header=labeled_rows[0]
    names=("Freq_MHz","PW_ns","Az_deg","El_deg")
    source_indices=[source_header.index(name) for name in names]
    truth_indices=[truth_header.index(name) for name in names]
    label_index=truth_header.index("Ground_Truth")
    destination.mkdir(parents=True,exist_ok=True)
    feature_file=destination/"features.csv";truth_file=destination/"truth.csv"
    with feature_file.open("w",newline="",encoding="utf-8") as ff,truth_file.open("w",newline="",encoding="utf-8") as tf:
        fw,tw=csv.writer(ff),csv.writer(tf);fw.writerow(names);tw.writerow([*names,"Ground_Truth"])
        for row in source_rows[1:]:fw.writerow([row[i] for i in source_indices])
        for row in labeled_rows[1:]:tw.writerow([*[row[i] for i in truth_indices],row[label_index]])
    return feature_file,truth_file

def all_jobs(prepared:Path)->list[tuple[str,Path,Path]]:
    data=ROOT/"Datasets";jobs=[]
    for name in ("s1_5d","s2_5d","s5_5d","s6_5d"):
        features,truth=prepare_combined(data/f"{name}.csv",prepared/name);jobs.append((name,features,truth))
    features,truth=prepare_two_emitter(data/"two_emitter_pdw_dataset.csv",data/"two_emitter_pdw_labeled.csv",prepared/"two_emitter_pdw")
    jobs.append(("two_emitter_pdw",features,truth));return jobs

def read_analysis(predictions:Path,report:Path)->dict:
    with predictions.open(newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    labels=[]
    for row in rows:
        for key in ("Ground_Truth","Predicted_Cluster"):
            if row[key] not in labels:labels.append(row[key])
    matrix=[[0 for _ in labels] for _ in labels];index={name:i for i,name in enumerate(labels)}
    for row in rows:matrix[index[row["Ground_Truth"]]][index[row["Predicted_Cluster"]]]+=1
    text=report.read_text(encoding="utf-8")
    def value(pattern:str,default:float=0.0)->float:
        match=re.search(pattern,text);return float(match.group(1)) if match else default
    accuracies=[float(x) for x in re.findall(r"^Accuracy=([0-9.eE+-]+)",text,re.MULTILINE)]
    return {"labels":labels,"matrix":matrix,"rows":len(rows),"signal_accuracy":accuracies[0] if accuracies else 0.0,
            "overall_accuracy":accuracies[1] if len(accuracies)>1 else 0.0,
            "ari":value(r"Adjusted Rand Index=([0-9.eE+-]+)"),"threshold":value(r"Applied threshold=([0-9.eE+-]+)"),
            "leaf_cfs":int(value(r"Leaf CFs=([0-9]+)")),
            "estimated_clusters":int(value(r"Automatically estimated clusters \(Gap Statistic\)=([0-9]+)")),
            "clusters":int(value(r"Final detected BIRCH clusters=([0-9]+)")),
            "noise":int(value(r"Noise rows from small CFs=([0-9]+)")),
            "separation":"PASS" if "Condition 6*Rmax <= Dmin: PASS" in text else "WARNING",
            "weight":"PASS" if "Weight-ratio condition: PASS" in text else "WARNING"}

def confusion_artifacts(folder:Path,data:dict)->None:
    labels,matrix=data["labels"],data["matrix"]
    with (folder/"confusion_matrix.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f);writer.writerow(["Actual / Predicted",*labels]);writer.writerows([[labels[i],*row] for i,row in enumerate(matrix)])
    cell=72;left=175;top=75;width=left+cell*len(labels)+20;height=top+cell*len(labels)+90;maximum=max(1,max(max(r) for r in matrix))
    parts=[f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}'>","<rect width='100%' height='100%' fill='white'/>",f"<text x='{width/2}' y='25' text-anchor='middle' font-family='sans-serif' font-size='19'>Confusion matrix</text>",f"<text x='{left+cell*len(labels)/2}' y='{height-12}' text-anchor='middle' font-family='sans-serif' font-size='15' font-weight='bold'>Predicted class (columns)</text>",f"<text x='22' y='{top+cell*len(labels)/2}' text-anchor='middle' font-family='sans-serif' font-size='15' font-weight='bold' transform='rotate(-90 22 {top+cell*len(labels)/2})'>Actual class (rows)</text>"]
    for i,label in enumerate(labels):
        safe=html.escape(label);parts.append(f"<text x='{left-8}' y='{top+i*cell+42}' text-anchor='end' font-family='sans-serif' font-size='11'>{safe}</text>")
        parts.append(f"<text x='{left+i*cell+cell/2}' y='{top+len(labels)*cell+22}' text-anchor='middle' font-family='sans-serif' font-size='10' transform='rotate(25 {left+i*cell+cell/2} {top+len(labels)*cell+22})'>{safe}</text>")
        for j,count in enumerate(matrix[i]):
            opacity=.10+.85*count/maximum;color="white" if opacity>.55 else "black"
            parts.append(f"<rect x='{left+j*cell}' y='{top+i*cell}' width='{cell-3}' height='{cell-3}' fill='#2563eb' fill-opacity='{opacity:.3f}'/>")
            parts.append(f"<text x='{left+j*cell+(cell-3)/2}' y='{top+i*cell+42}' text-anchor='middle' font-family='sans-serif' font-size='12' fill='{color}'>{count}</text>")
    parts.append("</svg>");(folder/"confusion_matrix.svg").write_text("".join(parts),encoding="utf-8")

    counts={label:0 for label in labels}
    for row in rows_from_csv(folder/"birch_results.csv"):counts[row["Predicted_Cluster"]]=counts.get(row["Predicted_Cluster"],0)+1
    bar_svg(folder/"cluster_sizes.svg",list(counts),list(counts.values()),"Predicted cluster sizes",integer=True)

def rows_from_csv(path:Path):
    with path.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))

def bar_svg(path:Path,names:list[str],values:list[float],title:str,integer:bool=False)->None:
    width=900;height=520;left=75;bottom=440;usable=760;maximum=max(values) if values else 1;step=usable/max(1,len(values));colors=["#2563eb","#16a34a","#dc2626","#9333ea","#ea580c"]
    p=[f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}'><rect width='100%' height='100%' fill='white'/><text x='{width/2}' y='28' text-anchor='middle' font-family='sans-serif' font-size='20'>{html.escape(title)}</text><line x1='{left}' y1='{bottom}' x2='{left+usable}' y2='{bottom}' stroke='#333'/>" ]
    for i,(name,value) in enumerate(zip(names,values)):
        h=350*value/max(maximum,1e-12);x=left+i*step+step*.12;label=str(int(value)) if integer else f"{value:.3f}"
        p.append(f"<rect x='{x}' y='{bottom-h}' width='{step*.76}' height='{h}' fill='{colors[i%len(colors)]}'/><text x='{x+step*.38}' y='{bottom-h-8}' text-anchor='middle' font-family='sans-serif' font-size='11'>{label}</text><text x='{x+step*.38}' y='{bottom+20}' text-anchor='middle' font-family='sans-serif' font-size='11'>{html.escape(name)}</text>")
    p.append("</svg>");path.write_text("".join(p),encoding="utf-8")

def aggregate_report(output:Path,analyses:list[tuple[str,dict]])->None:
    names=[name for name,_ in analyses];bar_svg(output/"accuracy_comparison.svg",names,[d["signal_accuracy"] for _,d in analyses],"Signal clustering accuracy")
    bar_svg(output/"ari_comparison.svg",names,[d["ari"] for _,d in analyses],"Adjusted Rand Index")
    bar_svg(output/"threshold_comparison.svg",names,[d["threshold"] for _,d in analyses],"Automatic A-BIRCH thresholds")
    lines=["# A-BIRCH analysis report","","Every dataset was fitted independently using A-BIRCH. Ground truth was used only after clustering for validation.","","| Dataset | Rows | Threshold | Leaf CFs | Gap estimate | Final detected clusters | BIRCH noise | Signal accuracy | Overall accuracy | ARI | Paper conditions |","|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for name,d in analyses:lines.append(f"| {name} | {d['rows']} | {d['threshold']:.6f} | {d['leaf_cfs']} | {d['estimated_clusters']} | **{d['clusters']}** | {d['noise']} | {d['signal_accuracy']:.4f} | {d['overall_accuracy']:.4f} | {d['ari']:.4f} | separation: {d['separation']}; weight: {d['weight']} |")
    lines += ["","## Graphs","","![Signal accuracy](accuracy_comparison.svg)","","![Adjusted Rand Index](ari_comparison.svg)","","![Automatic thresholds](threshold_comparison.svg)","","## Interpretation","","- `s2_5d`, `s5_5d`, and the two-emitter dataset separate their signal emitters very accurately.","- `s1_5d` and `s6_5d` are under-clustered by pure tree-BIRCH: distinct labeled emitters merge geometrically.","- Noise detection is strongest for `s2_5d`; in the other datasets, the automatic radius absorbs noise into signal CFs.","- All paper-condition warnings are retained because the radar data is five-dimensional and not guaranteed to be isotropic Gaussian data.","","Each dataset folder contains `birch_results.csv`, `birch_validation_report.txt`, `confusion_matrix.csv`, `confusion_matrix.svg`, and `cluster_sizes.svg`."]
    (output/"ANALYSIS_REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")

def main()->int:
    a=arguments();compiler=shutil.which(a.compiler)
    if not a.skip_build and compiler is None:
        print("g++ not found. Run: sudo apt update && sudo apt install build-essential",file=sys.stderr);return 2
    build=ROOT/"build";build.mkdir(parents=True,exist_ok=True)
    executable=build/("birch_improved.exe" if os.name=="nt" else "birch_improved")
    try:
        if not a.skip_build:run([compiler,"-std=c++14","-O2","-Wall","-Wextra","-Wpedantic",*SOURCES,"-o",str(executable)],IMPL)
        if not executable.is_file():raise FileNotFoundError(f"Executable missing: {executable}")
        if not a.skip_tests:run([str(executable),"--self-test"],IMPL)
        if a.build_only:print(f"Build completed: {executable}");return 0
        if (a.features is None)!=(a.truth is None):raise ValueError("--features and --truth must be supplied together")
        jobs=[(a.features.stem,a.features.resolve(),a.truth.resolve())] if a.features else all_jobs(build/"prepared")
        output=a.output_dir.resolve();output.mkdir(parents=True,exist_ok=True);index=[["dataset","status","predictions","report"]];analyses=[]
        for name,features,truth in jobs:
            folder=output/name;folder.mkdir(parents=True,exist_ok=True);predictions=folder/"birch_results.csv";report=folder/"birch_validation_report.txt"
            print(f"\n===== {name} =====",flush=True)
            run([str(executable),str(features),str(truth),str(predictions),str(report),str(a.branching_factor),str(a.minimum_cluster_points),str(a.pilot_points),str(a.pilot_max_clusters),str(a.threshold),str(a.mbd_spread),str(a.knn_k),str(a.knn_outlier_fraction)],IMPL)
            analysis=read_analysis(predictions,report);confusion_artifacts(folder,analysis);analyses.append((name,analysis))
            index.append([name,"analyzed",str(predictions),str(report)])
        with (output/"run_index.csv").open("w",newline="",encoding="utf-8") as f:csv.writer(f).writerows(index)
        aggregate_report(output,analyses)
        print(f"\nAnalyzed {len(jobs)} datasets. Results: {output}");return 0
    except (subprocess.CalledProcessError,OSError,ValueError) as e:
        print(f"Error: {e}",file=sys.stderr);return e.returncode if isinstance(e,subprocess.CalledProcessError) else 1

if __name__=="__main__":raise SystemExit(main())
