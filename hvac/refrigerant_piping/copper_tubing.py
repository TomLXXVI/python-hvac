from typing import cast
from dataclasses import dataclass
import shelve
from hvac import Quantity


Q_ = Quantity


@dataclass
class CopperTube:
    """Dataclass with specifications of copper tube referred to a nominal
    diameter.

    Parameters
    ----------
    DN:
        Nominal tube diameter.
    D_ext:
        Outside tube diameter.
    D_int:
        Inside tube diameter.
    t:
        Tube wall thickness.
    """
    DN: str
    D_ext: Quantity
    D_int: Quantity
    t: Quantity


class CopperTubing:
    """Class that saves and loads copper tube specifications to/from a shelf on
    disk.
    """
    db_path: str  # file path to the shelf (database, db)

    @staticmethod
    def _create_record(DN, D_ext, D_int, t):
        record = {
            'DN': DN,
            'D_ext': D_ext.to('inch').m,
            'D_int': D_int.to('inch').m,
            't': t.to('inch').m,
        }
        return record

    @staticmethod
    def _refurbish_record(record) -> CopperTube:
        copper_tube = CopperTube(
            DN=record['DN'],
            D_ext=Q_(record['D_ext'], 'inch'),
            D_int=Q_(record['D_int'], 'inch'),
            t=Q_(record['t'], 'inch'),
        )
        return copper_tube

    @classmethod
    def add_record(
        cls,
        DN: str,
        D_ext: Quantity,
        D_int: Quantity,
        t: Quantity
    ) -> None:
        """
        Adds a copper tube record to the copper tubing shelf.

        Parameters
        ----------
        DN:
            Nominal diameter of the tube, entered as a string. This will be
            used as the key to select a tube from the shelf.
        D_ext:
            Outside diameter.
        D_int:
            Inside diameter.
        t:
            Wall thickness.
        """
        with shelve.open(cls.db_path) as shelf:
            record = cls._create_record(DN, D_ext, D_int, t)
            shelf[DN] = record

    @classmethod
    def get_record(cls, DN: str) -> CopperTube | None:
        """Returns a single copper tube record from the copper tubing shelf."""
        with shelve.open(cls.db_path) as shelf:
            try:
                record = cast(dict, shelf[DN])
            except KeyError:
                return None
            return cls._refurbish_record(record)

    @classmethod
    def get_records(cls, *dns: str) -> tuple[CopperTube | None, ...]:
        """Get multiple copper tubes from the copper tubing shelf. When
        called without parameters, all copper tubes in the shelf are returned.
        """
        records = []
        with shelve.open(cls.db_path) as shelf:
            if dns:
                for dn in dns:
                    try:
                        record = cast(dict, shelf[dn])
                    except KeyError:
                        records.append(None)
                    records.append(cls._refurbish_record(record))
            else:
                records = [cls._refurbish_record(record) for record in shelf.values()]
        return tuple(records)
