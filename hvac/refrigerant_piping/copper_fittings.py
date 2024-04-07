import pandas as pd
from hvac import Quantity

Q_ = Quantity


class CopperFittings:
    db_path: str = ''
    table: pd.DataFrame = None

    @classmethod
    def _create_table(cls) -> None:
        """Creates a Pandas DataFrame with the equivalent lengths of fittings
        in feet. The index of the data frame points to nominal tube diameters;
        the columns list the different fittings or accessories:
        - 'long-radius-elbow',
        - 'short-radius-elbow',
        - 'tee-thru-flow', 'tee-branch-flow',
        - 'ball-valve', 'sight-glass'
        """
        DN = [
            '1/2', '5/8', '3/4', '7/8',
            '1 1/8', '1 3/8', '1 5/8',
            '2 1/8', '2 5/8', '3 1/8'
        ]
        fitting = [
            'long-radius-elbow', 'short-radius-elbow',
            'tee-thru-flow', 'tee-branch-flow',
            'ball-valve', 'sight-glass'
        ]
        L_eq = [
            [0.3, 0.4, 1.0, 1.0, 1.0, 1.0],
            [0.4, 0.5, 1.0, 1.3, 1.0, 1.0],
            [0.5, 0.6, 1.0, 1.4, 1.0, 1.0],
            [0.6, 0.7, 1.0, 1.8, 1.0, 1.0],
            [0.8, 0.9, 1.0, 2.6, 1.0, 1.0],
            [0.9, 1.4, 1.0, 3.4, 1.0, 1.0],
            [1.0, 1.7, 1.0, 4.4, 1.0, 1.0],
            [1.4, 2.3, 1.0, 5.9, 1.0, 1.0],
            [1.5, 3.4, 1.0, 8.0, 1.0, 1.0],
            [1.7, 5.0, 1.0, 10.1, 1.0, 1.0]
        ]
        df = pd.DataFrame(data=L_eq, index=DN, columns=fitting)
        cls._save_table(df)

    @classmethod
    def _save_table(cls, df):
        try:
            with open(cls.db_path, 'w') as fh:
                df.to_csv(fh)
        except FileNotFoundError:
            raise FileNotFoundError(
                'Unknown path: Please, set `db_path` '
                'to correct destination.'
            ) from None

    @classmethod
    def _load_table(cls):
        try:
            with open(cls.db_path) as fh:
                cls.table = pd.read_csv(fh, index_col=0)
        except FileNotFoundError:
            cls._create_table()
            cls._load_table()

    @classmethod
    def get_Leq(cls, fitting: str, DN: str) -> Quantity:
        """Returns the equivalent length of the indicated fitting with a nominal
        diameter `DN`.
        """
        if cls.table is None:
            cls._load_table()
        try:
            Leq = cls.table.loc[DN, fitting]
        except KeyError:
            raise KeyError(
                'Unknown value for `fitting` or `DN`: '
                'Check value for `fitting` or `DN`.'
            )
        return Q_(Leq, 'feet')
