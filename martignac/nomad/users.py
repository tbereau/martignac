import logging
from dataclasses import asdict, field

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

    def as_dict(self) -> dict:
        return asdict(self)


@ttl_cache(maxsize=128, ttl=180)
def get_user_by_id(user_id: str, use_prod: bool = True, timeout_in_sec: int = 10) -> NomadUser:
    logger.info(f"retrieving user {user_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(f"/users/{user_id}", timeout_in_sec=timeout_in_sec)
    user_schema = class_schema(NomadUser)
    return user_schema().load(response)
