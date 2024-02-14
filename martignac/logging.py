import logging

LOGGING_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

logging.basicConfig(format=LOGGING_FORMAT, datefmt="%y-%m-%d %H:%M:%S")
logger = logging.getLogger()
disable_modules = ["requests", "web3", "urllib3", "gql.transport", "parso"]
for module in disable_modules:
    logging.getLogger(module).setLevel(logging.WARNING)
