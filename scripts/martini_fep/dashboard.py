from signac_dashboard import Dashboard
from martignac.workflows.alchemical_transformation import AlchemicalTransformationFlow


class MyDashboard(Dashboard):
    pass


if __name__ == '__main__':
    MyDashboard(project=AlchemicalTransformationFlow()).main()
