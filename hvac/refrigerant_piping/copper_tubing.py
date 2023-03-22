from typing import cast
from dataclasses import dataclass
import shelve
from hvac import Quantity


Q_ = Quantity


@dataclass
class CopperTube:
    """Dataclass holding the specifications of a copper tube with given nominal
    diameter."""
    dn: str
    do: Quantity
    di: Quantity
    t: Quantity


class CopperTubing:
    """Class that saves and loads copper tube specifications to/from a shelf on
    disk."""
    db_path: str  # file path to the shelf (database, db)

    @staticmethod
    def _create_record(dn, do, di, t):
        record = {
            'dn': dn,
            'do': do.to('inch').m,
            'di': di.to('inch').m,
            't': t.to('inch').m,
        }
        return record

    @staticmethod
    def _refurbish_record(record):
        copper_tube = CopperTube(
            dn=record['dn'],
            do=Q_(record['do'], 'inch'),
            di=Q_(record['di'], 'inch'),
            t=Q_(record['t'], 'inch'),
        )
        return copper_tube

    @classmethod
    def add_record(cls, dn: str, do: Quantity, di: Quantity, t: Quantity) -> None:
        """
        Add copper tube record to copper tubing shelf.

        Parameters
        ----------
        dn:
            Nominal diameter of the tube, entered as a string. This will be
            used as the key to select a tube from the shelf.
        do:
            Outside diameter.
        di:
            Inside diameter.
        t:
            Wall thickness.
        """
        with shelve.open(cls.db_path) as shelf:
            record = cls._create_record(dn, do, di, t)
            shelf[dn] = record

    @classmethod
    def get_record(cls, dn: str) -> CopperTube | None:
        """Get single copper tube record from copper tubing shelf."""
        with shelve.open(cls.db_path) as shelf:
            try:
                record = cast(dict, shelf[dn])
            except KeyError:
                return None
            return cls._refurbish_record(record)

    @classmethod
    def get_records(cls, *dns: str) -> tuple[CopperTube | None, ...]:
        """Get multiple copper tube records from copper tubing shelf."""
        records = []
        with shelve.open(cls.db_path) as shelf:
            for dn in dns:
                try:
                    record = cast(dict, shelf[dn])
                except KeyError:
                    records.append(None)
                records.append(cls._refurbish_record(record))
        return tuple(records)
