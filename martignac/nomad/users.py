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
    logger.info(f"retrieving user {user_name} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(f"/users?prefix={user_name}", timeout_in_sec=timeout_in_sec).get("data", [])
    return [class_schema(NomadUser)().load(user) for user in response]


@ttl_cache(maxsize=128, ttl=180)
def get_user_by_id(user_id: str, use_prod: bool = True, timeout_in_sec: int = 10) -> NomadUser:
    logger.info(f"retrieving user {user_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(f"/users/{user_id}", timeout_in_sec=timeout_in_sec)
    user_schema = class_schema(NomadUser)
    return user_schema().load(response)


@ttl_cache(maxsize=128, ttl=180)
def who_am_i(use_prod: bool = True, timeout_in_sec: int = 10) -> NomadUser:
    logger.info(f"retrieving self user info on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request("/users/me", with_authentication=True, timeout_in_sec=timeout_in_sec)
    user_schema = class_schema(NomadUser)
    return user_schema().load(response)
