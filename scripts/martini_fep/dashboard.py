from signac_dashboard import Dashboard
from project import MartiniProject


class MyDashboard(Dashboard):
    pass


if __name__ == '__main__':
    MyDashboard(project=MartiniProject()).main()
