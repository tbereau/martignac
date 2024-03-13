import json
import logging
from typing import Any

import requests
from cachetools.func import ttl_cache
from decouple import config as environ

logger = logging.getLogger(__name__)

NOMAD_USERNAME = environ("NOMAD_USERNAME")
NOMAD_PASSWORD = environ("NOMAD_PASSWORD")
NOMAD_PROD_URL = "https://nomad-lab.eu/prod/v1/api/v1"
NOMAD_TEST_URL = "https://nomad-lab.eu/prod/v1/test/api/v1"
TIMEOUT_IN_SEC = 60


@ttl_cache(maxsize=128, ttl=180)
def get_authentication_token(
    use_prod: bool = False,
    username: str = NOMAD_USERNAME,
    password: str = NOMAD_PASSWORD,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
) -> str:
    url = NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL
    logger.info(f"Requesting authentication token @ {url}")
    response = requests.get(
        url + "/auth/token",
        params={"username": username, "password": password},
        timeout=timeout_in_sec,
    )
    if not response.status_code == 200:
        raise ValueError(f"Unexpected response {response.json()}")
    return response.json().get("access_token")


def get_nomad_request(
    section: str,
    use_prod: bool = False,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
    headers: dict = None,
    with_authentication: bool = False,
) -> json:
    url = NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL
    url += f"{'/' if section[0] != '/' else ''}{section}"
    logger.info(f"Sending get request @ {url}")
    if headers is None:
        headers = {}
    if with_authentication:
        token = get_authentication_token(use_prod=use_prod)
        headers |= {
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
        }
    response = requests.get(url, headers=headers, timeout=timeout_in_sec)
    if not response.status_code == 200:
        raise ValueError(f"Unexpected response {response.json()}")
    return response.json()


def get_nomad_base_url(use_prod: bool) -> str:
    return (NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL).removesuffix("/api/v1")


def post_nomad_request(
    section: str,
    headers: dict = None,
    data: Any = None,
    json_dict: dict = None,
    use_prod: bool = False,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
    with_authentication: bool = False,
) -> json:
    if headers is None:
        headers = {}
    if with_authentication:
        token = get_authentication_token(use_prod=use_prod)
        headers |= {
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
        }
    if data is None:
        data = {}
    if json_dict is None:
        json_dict = {}
    url = NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL
    url += f"{'/' if section[0] != '/' else ''}{section}"
    logger.info(f"Sending post request @ {url}")
    response = requests.post(url, headers=headers, json=json_dict, data=data, timeout=timeout_in_sec)
    if not response.status_code == 200:
        raise ValueError(f"Unexpected response {response.json()}")
    return response.json()


def delete_nomad_request(
    section: str,
    headers: dict = None,
    use_prod: bool = False,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
    with_authentication: bool = False,
) -> json:
    if headers is None:
        headers = {}
    if with_authentication:
        token = get_authentication_token(use_prod=use_prod)
        headers |= {
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
        }
    url = NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL
    url += f"{'/' if section[0] != '/' else ''}{section}"
    logger.info(f"Sending delete request @ {url}")
    response = requests.delete(url, headers=headers, timeout=timeout_in_sec)
    if not response.status_code == 200:
        raise ValueError(f"Unexpected response {response.json()}")
    return response.json()
