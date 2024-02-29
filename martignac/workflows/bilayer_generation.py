from flow import FlowProject

from martignac import config
from martignac.liquid_models.mixtures import LiquidMixture
from martignac.parsers.gromacs_topologies import Topology
from martignac.utils.gromacs import gromacs_simulation_command
from martignac.utils.insane import generate_bilayer_with_insane

conf = config()["bilayer_generation"]


class BilayerGenFlow(FlowProject):
    itp_files: list[str] = [str(e) for e in conf["itp_files"].get()]
    n_threads: int = conf["settings"]["n_threads"].get(int)
    box_length_xy: float = conf["settings"]["box_length_xy"].get(float)
    box_length_z: float = conf["settings"]["box_length_z"].get(float)
    solvent: LiquidMixture = LiquidMixture.from_list_of_dicts([s.get(dict) for s in conf["settings"]["solvent"]])
    minimize_mdp: str = conf["mdp_files"]["minimize"].get(str)
    equilibrate_mdp: str = conf["mdp_files"]["equilibrate"].get(str)
    system_name: str = conf["output_names"]["system"].get(str)
    gen_name: str = conf["output_names"]["states"]["generate"].get(str)
    min_name: str = conf["output_names"]["states"]["minimize"].get(str)
    equ_name: str = conf["output_names"]["states"]["equilibrate"].get(str)

    @classmethod
    def get_bilayer_gen_name(cls) -> str:
        return f"{cls.system_name}_{cls.gen_name}"

    @classmethod
    def get_bilayer_min_name(cls) -> str:
        return f"{cls.system_name}_{cls.min_name}"

    @classmethod
    def get_bilayer_equ_name(cls) -> str:
        return f"{cls.system_name}_{cls.equ_name}"

    @classmethod
    def get_bilayer_gen_gro(cls) -> str:
        return f"{cls.get_bilayer_gen_name()}.gro"

    @classmethod
    def get_bilayer_top(cls) -> str:
        return f"{cls.system_name}.top"

    @classmethod
    def get_bilayer_min_gro(cls) -> str:
        return f"{cls.get_bilayer_min_name()}.gro"

    @classmethod
    def get_bilayer_equ_gro(cls) -> str:
        return f"{cls.get_bilayer_equ_name()}.gro"


project = BilayerGenFlow().get_project()


@BilayerGenFlow.label
def generated_gen_gro(job):
    return job.isfile(BilayerGenFlow.get_bilayer_gen_gro())


@BilayerGenFlow.label
def generated_min_gro(job):
    return job.isfile(BilayerGenFlow.get_bilayer_min_gro())


@BilayerGenFlow.label
def generated_equ_gro(job):
    return job.isfile(BilayerGenFlow.get_bilayer_equ_gro())


@BilayerGenFlow.post(generated_gen_gro)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def generate_initial_bilayer(job):
    lipid_mixture = LiquidMixture.from_list_of_dicts(job.sp.lipids)
    return generate_bilayer_with_insane(
        lipids=lipid_mixture,
        solvent=BilayerGenFlow.solvent,
        box_length_xy=BilayerGenFlow.box_length_xy,
        box_length_z=BilayerGenFlow.box_length_z,
        gro_bilayer_gen=BilayerGenFlow.get_bilayer_gen_gro(),
        top_bilayer=BilayerGenFlow.get_bilayer_top(),
    )


@BilayerGenFlow.pre(generated_gen_gro)
@BilayerGenFlow.post(generated_min_gro)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def minimize(job):
    top = Topology.parse_top_file(BilayerGenFlow.get_bilayer_top())
    top.includes = BilayerGenFlow.itp_files
    top.output_top(BilayerGenFlow.get_bilayer_top())
    return gromacs_simulation_command(
        mdp=BilayerGenFlow.minimize_mdp,
        top=BilayerGenFlow.get_bilayer_top(),
        gro=BilayerGenFlow.get_bilayer_gen_gro(),
        name=BilayerGenFlow.get_bilayer_min_name(),
        n_threads=BilayerGenFlow.n_threads,
    )


@BilayerGenFlow.pre(generated_min_gro)
@BilayerGenFlow.post(generated_equ_gro)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    return gromacs_simulation_command(
        mdp=BilayerGenFlow.equilibrate_mdp,
        top=BilayerGenFlow.get_bilayer_top(),
        gro=BilayerGenFlow.get_bilayer_min_gro(),
        name=BilayerGenFlow.get_bilayer_equ_name(),
        n_threads=BilayerGenFlow.n_threads,
    )
