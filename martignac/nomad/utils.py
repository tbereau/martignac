import os
import json
import logging
import time
from functools import wraps
from typing import Any, Optional

import requests
from cachetools.func import ttl_cache
from decouple import AutoConfig
from decouple import config as environ
from requests_toolbelt.multipart.encoder import MultipartEncoder

logger = logging.getLogger(__name__)

environ = AutoConfig(search_path=os.environ["MARTIGNACDIR"])
NOMAD_USERNAME = environ("NOMAD_USERNAME")
NOMAD_PASSWORD = environ("NOMAD_PASSWORD")
NOMAD_PROD_URL = "https://nomad-lab.eu/prod/v1/api/v1"
NOMAD_TEST_URL = "https://nomad-lab.eu/prod/v1/test/api/v1"
TIMEOUT_IN_SEC = 60
NOMAD_SLEEP_INTERVAL_IN_SECONDS = 0.2


def rate_limiter(min_interval):
    """
    A decorator to enforce a minimum time interval between calls to the decorated function.

    Parameters:
        min_interval (float): Minimum time interval between consecutive calls in seconds.
    """

    def decorator(func):
        last_call = [0]  # Use a mutable object to keep track of the last call time

        @wraps(func)
        def wrapper(*args, **kwargs):
            elapsed = time.time() - last_call[0]
            if elapsed < min_interval:
                time.sleep(min_interval - elapsed)
            result = func(*args, **kwargs)
            last_call[0] = time.time()
            return result

        return wrapper

    return decorator


@ttl_cache(maxsize=128, ttl=180)
def get_authentication_token(
    use_prod: bool = False,
    username: str = NOMAD_USERNAME,
    password: str = NOMAD_PASSWORD,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
) -> str:
    """
    Retrieves an authentication token from the NOMAD API.

    This function requests an authentication token by providing user credentials. It supports both production
    and test environments of the NOMAD API. The function raises an exception if the request fails or if the
    response status code is not 200.

    Parameters:
        use_prod (bool): Flag to determine whether to use the production URL or the test URL. Defaults to False.
        username (str): The username for authentication. Defaults to NOMAD_USERNAME from environment variables.
        password (str): The password for authentication. Defaults to NOMAD_PASSWORD from environment variables.
        timeout_in_sec (int): The timeout for the request in seconds. Defaults to TIMEOUT_IN_SEC from environment variables.

    Returns:
        str: The authentication token as a string.

    Raises:
        ValueError: If the response from the server is not successful (status code 200).
    """
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


@rate_limiter(min_interval=NOMAD_SLEEP_INTERVAL_IN_SECONDS)
def get_nomad_request(
    section: str,
    use_prod: bool = False,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
    headers: Optional[dict] = None,
    with_authentication: bool = False,
    return_json: bool = True,
    accept_field: str = "application/json",
) -> Any:
    """
    Sends a GET request to the NOMAD API for a specified section of the database.

    This function constructs and sends a GET request to either the production or test environment of the NOMAD API,
    based on the provided parameters. It can optionally include authentication in the request headers. The response
    can be returned either as JSON or raw content, based on the `return_json` parameter.

    Parameters:
        section (str): The specific section of the NOMAD API to query.
        use_prod (bool): Determines whether to use the production or test URL for the NOMAD API. Defaults to False.
        timeout_in_sec (int): The timeout for the request in seconds. Defaults to TIMEOUT_IN_SEC from environment variables.
        headers (dict, optional): Additional headers to include in the request. Defaults to None.
        with_authentication (bool): If True, includes an authentication token in the request headers. Defaults to False.
        return_json (bool): If True, returns the response in JSON format; otherwise, returns raw content. Defaults to True.
        accept_field (str): The value for the 'Accept' header in the request, indicating the desired response format. Defaults to "application/json".

    Returns:
        Any: The response from the NOMAD API, formatted as JSON or raw content based on the `return_json` parameter.

    Raises:
        ValueError: If the response from the NOMAD API is not successful (status code 200).
    """
    url = NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL
    url += f"{'/' if section[0] != '/' else ''}{section}"
    logger.info(f"Sending get request @ {url}")
    if headers is None:
        headers = {}
    if with_authentication:
        token = get_authentication_token(use_prod=use_prod)
        headers |= {
            "Authorization": f"Bearer {token}",
            "Accept": accept_field,
        }
    response = requests.get(url, headers=headers, timeout=timeout_in_sec)
    if not response.status_code == 200:
        raise ValueError(f"Unexpected response {response.json()}")
    if return_json:
        return response.json()
    return response.content


def get_nomad_base_url(use_prod: bool) -> str:
    """
    Returns the base URL for the NOMAD API depending on the environment.

    This function provides the base URL for either the production or test environment of the NOMAD API. It is useful
    for constructing full API request URLs based on the desired environment.

    Parameters:
        use_prod (bool): A flag indicating whether to return the production URL. If False, the test URL is returned.

    Returns:
        str: The base URL for the NOMAD API.
    """
    return (NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL).removesuffix("/api/v1")


@rate_limiter(min_interval=NOMAD_SLEEP_INTERVAL_IN_SECONDS)
def post_nomad_request(
    section: str,
    headers: Optional[dict] = None,
    data: Any = None,
    json_dict: Optional[dict] = None,
    use_prod: bool = False,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
    with_authentication: bool = False,
) -> json:
    """
    Sends a POST request to the NOMAD API for a specified section of the database.

    This function constructs and sends a POST request to either the production or test environment of the NOMAD API,
    based on the provided parameters. It can optionally include authentication in the request headers and allows for
    sending data either as form data or JSON. The response is expected to be in JSON format.

    Parameters:
        section (str): The specific section of the NOMAD API to target with the POST request.
        headers (dict, optional): Additional headers to include in the request. Defaults to None.
        data (Any, optional): Form data to send with the request. Defaults to None.
        json_dict (dict, optional): JSON data to send with the request. Defaults to None.
        use_prod (bool): Determines whether to use the production or test URL for the NOMAD API. Defaults to False.
        timeout_in_sec (int): The timeout for the request in seconds. Defaults to TIMEOUT_IN_SEC from environment variables.
        with_authentication (bool): If True, includes an authentication token in the request headers. Defaults to False.

    Returns:
        json: The response from the NOMAD API in JSON format.

    Raises:
        ValueError: If the response from the NOMAD API is not successful (status code 200).
    """
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
    response = requests.post(
        url, headers=headers, json=json_dict, data=data, timeout=timeout_in_sec
    )
    if not response.status_code == 200:
        raise ValueError(f"Unexpected response {response.json()}")
    return response.json()


@rate_limiter(min_interval=NOMAD_SLEEP_INTERVAL_IN_SECONDS)
def delete_nomad_request(
    section: str,
    headers: Optional[dict] = None,
    use_prod: bool = False,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
    with_authentication: bool = False,
) -> json:
    """
    Sends a DELETE request to the NOMAD API for a specified section of the database.

    This function constructs and sends a DELETE request to either the production or test environment of the NOMAD API,
    based on the provided parameters. It can optionally include authentication in the request headers. The response is
    expected to be in JSON format.

    Parameters:
        section (str): The specific section of the NOMAD API to target with the DELETE request.
        headers (dict, optional): Additional headers to include in the request. Defaults to None.
        use_prod (bool): Determines whether to use the production or test URL for the NOMAD API. Defaults to False.
        timeout_in_sec (int): The timeout for the request in seconds. Defaults to TIMEOUT_IN_SEC from environment variables.
        with_authentication (bool): If True, includes an authentication token in the request headers. Defaults to False.

    Returns:
        json: The response from the NOMAD API in JSON format.

    Raises:
        ValueError: If the response from the NOMAD API is not successful (status code 200).
    """
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


@rate_limiter(min_interval=NOMAD_SLEEP_INTERVAL_IN_SECONDS)
def put_nomad_request(
    section: str,
    file: str,
    remote_path_to_file: str,
    headers: Optional[dict] = None,
    use_prod: bool = False,
    timeout_in_sec: int = TIMEOUT_IN_SEC,
    with_authentication: bool = False,
):
    if headers is None:
        headers = {}
    with open(file, "rb") as f:
        mp_encoder = MultipartEncoder(
            fields={"file": (remote_path_to_file, f, "application/json")}
        )

    if with_authentication:
        token = get_authentication_token(use_prod=use_prod)
        headers |= {
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
            "Content-Type": mp_encoder.content_type,
        }
    url = NOMAD_PROD_URL if use_prod else NOMAD_TEST_URL
    url += f"{'/' if section[0] != '/' else ''}{section}"
    logger.info(f"Sending put request @ {url}")
    response = requests.put(
        url, headers=headers, data=mp_encoder, timeout=timeout_in_sec
    )
    if not response.status_code == 200:
        raise ValueError(f"Unexpected response {response.json()}")
    return response.json()
