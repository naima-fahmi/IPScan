import os, time, csv
import configparser
import logging
import datetime
from methods.print_and_log import print_and_log
from methods.methods import make_full_dictionary
from methods.methods import generate_ipscan_events


# Read in the configuration values
config = configparser.ConfigParser()
config.read("./config/config.ini")

# Parameters
species = config["PARAMETERS"]["species"]
annotation_file = config["PARAMETERS"]["annotation_file"]
mode = config["PARAMETERS"]["mode"]
cores = config["PARAMETERS"]["cores"]
sample1 = config["PARAMETERS"]["sample1"]
sample2 = config["PARAMETERS"]["sample2"]
output_dir = config["PARAMETERS"]["output_dir"]

# Logging Folder
log_dir = config["LOGGING"]["log_dir"]

# Configure the logging object
log_file_name = os.path.join(log_dir, f"ipscan.{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
logging.basicConfig(filename=log_file_name, level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger()


# Begin the pipeline
print_and_log("----------------------------------------------------------------", logger)
print_and_log("| Beginning IPScan Run...                               |", logger)
print_and_log("----------------------------------------------------------------", logger)
logger.handlers[0].stream.write("\n\n")
print()

# Print the config file for logging
print_and_log("----------------------------------------------------------------", logger)
print_and_log("| Configuration File:                                          |", logger)
print_and_log("----------------------------------------------------------------", logger)
for section in config:
    if section == "DEFAULT":
        continue
    print_and_log(section, logger)
    for key in config[section]:
        print_and_log(f"{key}: {config[section][key]}", logger)
    print()
    logger.handlers[0].stream.write("\n")
print_and_log("----------------------------------------------------------------", logger)
logger.handlers[0].stream.write("\n\n")
print()


# Track time
start = time.perf_counter()

# generate IPA events list
prestart = time.perf_counter()
generate_ipscan_events(species, annotation_file, sample1, output_dir, logger)

"""
df = pd.read_csv(outdir+"IPScan_"+sample+".csv", delimiter = '\t')
df.to_excel(outdir+"IPScan_"+sample+".xlsx", sheet_name='Sheet1', index = None, header=True)
"""
print_and_log(f"Step 0 Completed in {time.perf_counter() - prestart:.2f} seconds", logger)
print_and_log("----------------------------------------------------------------", logger)
logger.handlers[0].stream.write("\n\n")
print()


# Begin the pipeline
print_and_log("----------------------------------------------------------------", logger)
print_and_log(f"| IPScan pipeline completed in {time.perf_counter() - start:.2f} seconds                 |", logger)
print_and_log("----------------------------------------------------------------", logger)

