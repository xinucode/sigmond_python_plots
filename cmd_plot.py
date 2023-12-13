import argparse

def pickup_config():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file")
    args = parser.parse_args()
    config_file = args.config
    return config_file
    
def pickup_config_and_info():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file")
    parser.add_argument("data", help="data file")
    args = parser.parse_args()
    config_file = args.config
    data_file = args.data
    return config_file, data_file