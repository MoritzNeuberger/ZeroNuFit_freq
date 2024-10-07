    import json
import numpy as np
from legendmeta import LegendMetadata, to_datetime
import argparse
from datetime import datetime, timezone


lmeta = LegendMetadata("/mnt/atlas01/users/neuberger/L200_stat_analysis/sandbox/legend-metadata/")
chmap = lmeta.channelmap()

def recursive_items(dictionary):
    """
    Function to iterate over a nested dict returning the key and value tuple for each dict and subdict.  
    """
    for key, value in dictionary.items():
        if type(value) is dict:
            yield (key, value)
            yield from recursive_items(value)
        else:
            yield (key, value)
            
def define_global_partition_spans(l_input):
    """
    Generates list of tuples defining the start and end timestamp of global paritions.
    """
    output = list() 
    for key, value in recursive_items(l_input):
        if key == "span_in_utc_s":
            output.append(value)
    time_stamps = np.sort(np.unique(np.concatenate(output)))
    return np.array([[time_stamps[i],time_stamps[i+1]] for i in range(len(time_stamps) - 1)])

def get_partition_name(parition_data,span):
    """
    Returns the name (key) of the local partition of a specific detector falling into the given time span.
    """
    for key, value in parition_data.items():
        if key == "default":
            continue
        if not "span_in_utc_s" in value.keys():
            continue
        if value["span_in_utc_s"][0] <= span[0] and value["span_in_utc_s"][1] >= span[1]:
            return key
    return None

def get_run_range(span,buffer=10):
    """
    Returns a list of dicts containing period and run information falling into the given time span.
    """
    runs = []
    for period in lmeta.dataprod.config.analysis_runs.keys():
        for run in lmeta.dataprod.config.analysis_runs[period]:
            if "phy" in lmeta.dataprod.runinfo[period][run].keys():
                start_key_ts = int(to_datetime(lmeta.dataprod.runinfo[period][run]["phy"]["start_key"]).replace(tzinfo=timezone.utc).timestamp()) 
                if (start_key_ts >= span[0] - buffer) & (start_key_ts < span[1] - buffer):
                    runs.append({"period":period,"run":run})
    return runs

def get_exposure_of_det_in_span(det_name,span):
    """
    Returns the exposure of a detector in a specific time span.
    """
    with open("skip_runs_of_dets.json","r") as f:
        skip_runs_of_dets = json.load(f)
    skip_runs = []
    if det_name in skip_runs_of_dets.keys():
        skip_runs = skip_runs_of_dets[det_name]
    runs_in_span = get_run_range(span)
    if len(runs_in_span) == 0:
        return 0
    livetime = 0
    for run_info in runs_in_span:
        skip_flag = False
        for skip in skip_runs:
            #print(det_name,run_info["period"],skip["period"],run_info["run"],skip["run"],run_info["period"] == skip["period"] , run_info["run"] == skip["run"])
            if run_info["period"] == skip["period"] and run_info["run"] == skip["run"]:
                skip_flag = True
        if skip_flag:
            continue
        chmap_local = lmeta.channelmap(lmeta.dataprod.runinfo[run_info["period"]][run_info["run"]]["phy"]["start_key"])

        livetime += lmeta.dataprod.runinfo[run_info["period"]][run_info["run"]]["phy"]["livetime_in_s"] * (chmap_local[det_name]["analysis"]["usability"] == "on")
    return livetime / (3600 * 24 * 365.25) * chmap[det_name]["production"]["mass_in_g"] * 1e-3

def calculate_total_efficiency(global_default,detector_default,partition_specific,det_name):
    """
    Calculates the total efficiency and uncertainty.
    """
    terms = [
                partition_specific["ovbb_acceptance"]["psd"],
                partition_specific["ovbb_acceptance"]["quality"],
                detector_default["ovbb_acceptance"]["active_volume"],
                detector_default["ovbb_acceptance"]["containment"],
                global_default["ovbb_acceptance"]["lar"],
                lmeta.hardware.detectors.germanium.diodes[det_name].production.enrichment
            ]
    eff = np.prod([entry["val"] for entry in terms])
    eff_err = np.sqrt( np.sum([ (terms[i]["unc"] * np.prod([ terms[j]["val"] for j in range(len(terms)) if i != j ]))**2 for i in range(len(terms)) ]) )
    return eff, eff_err

def generate_one_partition_data(l_input,span):
    """
    Generates the information for one global partition.
    """
    output = {}
    for key, value in l_input.items():
        if key == "default":
            continue
        partition_name = get_partition_name(value,span)
        if not partition_name:
            continue
        exposure = get_exposure_of_det_in_span(key,span)
        if exposure == 0:
            continue
        eff, eff_err = calculate_total_efficiency(l_input["default"]["default"],
                                                  value["default"],
                                                  value[partition_name],
                                                 key)
        output[key] = {
            "detnum": chmap[key].daq.rawid,
            "detname": key,
            "fwhm": value[partition_name]["fwhm_in_keV"]["val"],
            "fwhm_err": value[partition_name]["fwhm_in_keV"]["unc"],
            "eff": eff,
            "eff_err": eff_err,
            "exposure": exposure,
            "shift": value[partition_name]["energy_bias_in_keV"]["val"],
            "shift_err": value[partition_name]["energy_bias_in_keV"]["unc"],
            "bkgname": "l200_nu24_BI",
            "sysnam": "l200_sys1" 
        }
    if len(output.keys()) > 0:
        return output
    else:
        return None

def generate_global_partition_data(l_input):
    """
    Generates the information for all global paritions.
    """
    spans = define_global_partition_spans(l_input)
    out_spans = []
    span_data = []
    for span in spans:
        tmp = generate_one_partition_data(l_input,span)
        if tmp:
            out_spans.append(span)
            span_data.append(generate_one_partition_data(l_input,span))
    return {"spans": out_spans, "span_data":span_data}

def write_output(global_partition_data,file_name,iter_offset):
    """
    Writes the global partition information into a text file in GERDA style.
    """
    
    def generate_det_number_list(global_partition_data,iter_offset):
        det_number_dict = {}
        iterator = iter_offset
        for part in global_partition_data["span_data"]:
            for key in part.keys():
                if not key in det_number_dict:
                    det_number_dict[key] = iterator
                    iterator+=1
        return det_number_dict
    
    
    det_num_dict = generate_det_number_list(global_partition_data,iter_offset)
    
    with open(file_name,"w") as f:
        for i in range(len(global_partition_data["spans"])):
            f.write("t {:d} {:d}\n".format(global_partition_data["spans"][i][0],global_partition_data["spans"][i][1]))
            for det_name in global_partition_data["span_data"][i].keys():
                tmp_data = global_partition_data["span_data"][i][det_name]
                f.write("{:d} {} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {} {}\n".format(
                    det_num_dict[tmp_data["detname"]],
                    tmp_data["detname"],
                    tmp_data["fwhm"],
                    tmp_data["fwhm_err"],
                    tmp_data["eff"],
                    tmp_data["eff_err"],
                    tmp_data["exposure"],
                    tmp_data["shift"],
                    tmp_data["shift_err"],
                    tmp_data["bkgname"],
                    tmp_data["sysnam"]
                ))
                
parser = argparse.ArgumentParser(description='Converts a parameter file for statistical analysis from LEGEND to GERDA style.')
parser.add_argument('input', help='LEGEND style parameter file')
parser.add_argument('output', help='GERDA style parameter file')
parser.add_argument('--iter_offset',help="Detector iterator offset",default=0)
args = vars(parser.parse_args())
                
with open(args["input"],"r") as f:
    l_input = json.load(f)  
    
global_partition_data = generate_global_partition_data(l_input)
write_output(global_partition_data,args["output"],iter_offset=args["iter_offset"])
