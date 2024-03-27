import datetime as dt
import logging
from typing import Any, Optional

from cachetools.func import ttl_cache
from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass

from martignac.nomad.users import NomadUser, get_user_by_id
from martignac.nomad.utils import delete_nomad_request, get_nomad_base_url, get_nomad_request, post_nomad_request

logger = logging.getLogger(__name__)


class NomadUploadSchema(Schema):
    @pre_load
    def convert_users(self, data, **kwargs):
        data["main_author"] = get_user_by_id(user_id=data["main_author"]).as_dict()
        data["writers"] = [get_user_by_id(user_id=w).as_dict() for w in data["writers"]]
        data["reviewers"] = [get_user_by_id(user_id=r).as_dict() for r in data["reviewers"]]
        data["viewers"] = [get_user_by_id(user_id=v).as_dict() for v in data["viewers"]]
        return data


@dataclass(frozen=True)
class NomadUpload:
    upload_id: str
    upload_create_time: dt.datetime
    main_author: NomadUser
    process_running: bool
    current_process: str
    process_status: str
    last_status_message: str
    errors: list[Any]
    warnings: list[Any]
    coauthors: list[str]
    coauthor_groups: list[Any]
    reviewers: list[NomadUser]
    reviewer_groups: list[Any]
    writers: list[NomadUser]
    writer_groups: list[Any]
    viewers: list[NomadUser]
    viewer_groups: list[Any]
    published: bool
    published_to: list[Any]
    with_embargo: bool
    embargo_length: float
    license: str
    entries: int
    n_entries: Optional[int] = None
    upload_files_server_path: Optional[str] = None
    publish_time: Optional[dt.datetime] = None
    references: Optional[list[str]] = None
    datasets: Optional[list[str]] = None
    external_db: Optional[str] = None
    upload_name: Optional[str] = None
    comment: Optional[str] = None
    use_prod: Optional[bool] = None
    complete_time: Optional[dt.datetime] = None

    @property
    def base_url(self) -> Optional[str]:
        if self.use_prod is not None:
            return get_nomad_base_url(self.use_prod)
        return None

    @property
    def nomad_gui_url(self) -> str:
        if self.base_url is None:
            raise ValueError(f"missing attribute 'use_prod' for upload {self}")
        return f"{self.base_url}/gui/user/uploads/upload/id/{self.upload_id}"


@ttl_cache(maxsize=128, ttl=180)
def get_all_my_uploads(use_prod: bool = False, timeout_in_sec: int = 10) -> list[NomadUpload]:
    logger.info(f"retrieving all uploads on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        "/uploads",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return [upload_class_schema().load({**r, "use_prod": use_prod}) for r in response["data"]]


def get_upload_by_id(upload_id: str, use_prod: bool = False, timeout_in_sec: int = 10) -> NomadUpload:
    logger.info(f"retrieving upload {upload_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        f"/uploads/{upload_id}",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return upload_class_schema().load({**response["data"], "use_prod": use_prod})


def delete_upload(upload_id: str, use_prod: bool = False, timeout_in_sec: int = 10) -> NomadUpload:
    logger.info(f"deleting upload {upload_id} on {'prod' if use_prod else 'test'} server")
    response = delete_nomad_request(
        f"/uploads/{upload_id}",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return upload_class_schema().load({**response["data"], "use_prod": use_prod})


def upload_files_to_nomad(filename: str, use_prod: bool = False, timeout_in_sec: int = 30) -> str:
    logger.info(f"uploading file {filename} on {'prod' if use_prod else 'test'} server")
    with open(filename, "rb") as f:
        response = post_nomad_request(
            "/uploads",
            with_authentication=True,
            data=f,
            use_prod=use_prod,
            timeout_in_sec=timeout_in_sec,
        )
    upload_id = response.get("upload_id")
    if upload_id:
        logger.info(f"successful upload to {upload_id}")
        return upload_id
    else:
        logger.error(f"could not upload {filename}. Response {response}")


def publish_upload(upload_id: str, use_prod: bool = False, timeout_in_sec: int = 10) -> dict:
    logger.info(f"publishing upload {upload_id} on {'prod' if use_prod else 'test'} server")
    response = post_nomad_request(
        f"/uploads/{upload_id}/action/publish",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    return response


def edit_upload_metadata(
    upload_id: str,
    upload_name: Optional[str] = None,
    references: Optional[list[str]] = None,
    dataset_id: Optional[str] = None,
    embargo_length: Optional[float] = None,
    coauthors_ids: Optional[list[str]] = None,
    comment: Optional[str] = None,
    use_prod: bool = False,
    timeout_in_sec: int = 10,
) -> dict:
    logger.info(f"editing the metadata for upload {upload_id} on {'prod' if use_prod else 'test'} server")
    metadata = {"metadata": {}}
    if upload_name:
        metadata["metadata"]["upload_name"] = upload_name
    if references:
        metadata["metadata"]["references"] = references
    if dataset_id:
        metadata["metadata"]["datasets"] = dataset_id
    if embargo_length:
        metadata["metadata"]["embargo_length"] = embargo_length
    if coauthors_ids:
        metadata["metadata"]["coauthors"] = coauthors_ids
    if comment:
        metadata["metadata"]["comment"] = comment
    response = post_nomad_request(
        f"/uploads/{upload_id}/edit",
        with_authentication=True,
        json_dict=metadata,
        timeout_in_sec=timeout_in_sec,
    )
    return response
