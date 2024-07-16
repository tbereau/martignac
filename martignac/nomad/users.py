import datetime as dt
import logging
from dataclasses import asdict, field
from typing import Optional

from cachetools.func import ttl_cache
from marshmallow_dataclass import class_schema, dataclass

from martignac.nomad.utils import get_nomad_request

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class NomadUser:
    """
    Represents a user in the NOMAD system.

    This class defines the structure for storing user information as retrieved from the NOMAD system. It includes
    both mandatory and optional fields, with some fields being hidden from the representation for privacy or
    security reasons.

    Attributes:
        user_id (str): A unique identifier for the user, not shown in the representation.
        name (str): The full name of the user.
        first_name (str): The user's first name, not shown in the representation.
        last_name (str): The user's last name, not shown in the representation.
        username (str): The user's username, not shown in the representation.
        affiliation (str): The user's affiliated institution or organization, not shown in the representation.
        affiliation_address (str): The address of the affiliated institution or organization, not shown in the representation.
        email (Optional[str]): The user's email address, optional and not shown in the representation.
        is_oasis_admin (Optional[bool]): Flag indicating if the user is an OASIS admin, optional and not shown in the representation.
        is_admin (Optional[bool]): Flag indicating if the user is an admin, optional and not shown in the representation.
        repo_user_id (Optional[str]): A repository-specific user identifier, optional and not shown in the representation.
        created (Optional[dt.datetime]): The date and time when the user was created, optional and not shown in the representation.
    """

    user_id: str = field(repr=False)
    name: str
    first_name: str = field(repr=False)
    last_name: str = field(repr=False)
    username: str = field(repr=False)
    affiliation: str = field(repr=False)
    affiliation_address: str = field(repr=False)
    email: Optional[str] = field(repr=False, default=None)
    is_oasis_admin: Optional[bool] = field(repr=False, default=None)
    is_admin: Optional[bool] = field(repr=False, default=None)
    repo_user_id: Optional[str] = field(repr=False, default=None)
    created: Optional[dt.datetime] = field(repr=False, default=None)

    def as_dict(self) -> dict:
        return asdict(self)


@ttl_cache(maxsize=128, ttl=180)
def search_users_by_name(user_name: str, use_prod: bool = True, timeout_in_sec: int = 10) -> NomadUser:
    """
    Searches for users in the NOMAD system by name.

    This function queries the NOMAD system for users whose names match the provided prefix. It leverages caching to
    reduce the number of requests made for the same query. The search can be directed to either the production or
    test environment of the NOMAD system. A timeout can be specified to limit the duration of the request.

    Args:
        user_name (str): The prefix of the user name to search for.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to True.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        NomadUser: An instance or list of `NomadUser` objects matching the search criteria.

    Note:
        The function is decorated with `@ttl_cache` to cache the results for 180 seconds and limit the cache size to 128 entries.
    """
    logger.info(f"retrieving user {user_name} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(f"/users?prefix={user_name}", timeout_in_sec=timeout_in_sec).get("data", [])
    return [class_schema(NomadUser)().load(user) for user in response]


@ttl_cache(maxsize=128, ttl=180)
def get_user_by_id(user_id: str, use_prod: bool = True, timeout_in_sec: int = 10) -> NomadUser:
    """
    Retrieves a user's information from the NOMAD system by their unique user ID.

    This function makes a request to the NOMAD system to fetch the details of a user specified by their unique identifier.
    It allows the choice between querying the production or test environment of the NOMAD system. The function also
    supports specifying a timeout for the request. The results are cached to optimize performance and reduce the load
    on the NOMAD system.

    Args:
        user_id (str): The unique identifier of the user to retrieve.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to True.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        NomadUser: An instance of `NomadUser` containing the retrieved user's information.

    Note:
        The function is decorated with `@ttl_cache` to cache the results for 180 seconds and limit the cache size to 128 entries.
    """
    logger.info(f"retrieving user {user_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(f"/users/{user_id}", timeout_in_sec=timeout_in_sec)
    user_schema = class_schema(NomadUser)
    return user_schema().load(response)


@ttl_cache(maxsize=128, ttl=180)
def who_am_i(use_prod: bool = True, timeout_in_sec: int = 10) -> NomadUser:
    """
    Retrieves the information of the currently authenticated user from the NOMAD system.

    This function makes a request to the NOMAD system to fetch the details of the currently authenticated user.
    It allows specifying whether to query the production or test environment of the NOMAD system. Additionally,
    a timeout for the request can be defined to manage how long the function waits for a response from the NOMAD system.
    The results are cached to improve performance and reduce the load on the NOMAD system.

    Args:
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to True.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        NomadUser: An instance of `NomadUser` containing the information of the currently authenticated user.

    Note:
        The function is decorated with `@ttl_cache` to cache the results for 180 seconds and limit the cache size to 128 entries.
    """
    logger.info(f"retrieving self user info on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request("/users/me", with_authentication=True, timeout_in_sec=timeout_in_sec)
    user_schema = class_schema(NomadUser)
    return user_schema().load(response)
