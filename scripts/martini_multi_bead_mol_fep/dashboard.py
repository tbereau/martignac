from signac_dashboard import Dashboard

from martignac.workflows.free_energy_calculations.alchemical_transformation import AlchemicalTransformationFlow


class MyDashboard(Dashboard):
    pass


if __name__ == "__main__":
    MyDashboard(project=AlchemicalTransformationFlow()).main()
