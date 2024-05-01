from signac_dashboard import Dashboard

from martignac.workflows.system_generation.solute_solvation import SoluteSolvationFlow


class MyDashboard(Dashboard):
    pass


if __name__ == "__main__":
    MyDashboard(project=SoluteSolvationFlow()).main()
