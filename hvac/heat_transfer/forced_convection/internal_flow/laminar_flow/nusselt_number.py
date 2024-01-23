"""
Implementation of the correlations for fully developed and average Nusselt number
in case of laminar flow.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
from scipy.interpolate import interp1d
from hvac import Quantity


class CircularTube:

    class UniformHeatFlux:

        @staticmethod
        def fully_developed_local_nusselt_number() -> float:
            return 4.36

        @staticmethod
        def average_nusselt_number(Pr: float, Gz: float) -> float:
            """Average Nusselt number for simultaneously developing flow."""
            UHF = CircularTube.UniformHeatFlux
            Nu_fd = UHF.fully_developed_local_nusselt_number()
            Nu_avg = Nu_fd
            n = (0.1156 + 0.08569 / (Pr ** 0.4)) * Gz
            d = 1 + 0.1158 * (Gz ** 0.6)
            Nu_avg += n / d
            return Nu_avg

    class UniformWallTemperature:

        @staticmethod
        def fully_developed_local_nusselt_number() -> float:
            return 3.66

        @staticmethod
        def average_nusselt_number(Pr: float, Gz: float) -> float:
            """Average Nusselt number for simultaneously developing flow."""
            UWT = CircularTube.UniformWallTemperature
            Nu_fd = UWT.fully_developed_local_nusselt_number()
            Nu_avg = Nu_fd
            n = (0.049 + 0.020 / Pr) * (Gz ** 1.12)
            d = 1 + 0.065 * (Gz ** 0.7)
            Nu_avg += n / d
            return Nu_avg


class RectangularTube:

    class UniformHeatFlux:

        @staticmethod
        def fully_developed_local_nusselt_number(a: Quantity, b: Quantity) -> float:
            a = a.to('m').m
            b = b.to('m').m
            ar = min(a, b) / max(a, b)  # aspect ratio
            Nu_fd = 8.235 * (
                1
                - 2.042 * ar
                + 3.085 * (ar ** 2)
                - 2.477 * (ar ** 3)
                + 1.058 * (ar ** 4)
                - 0.186 * (ar ** 5)
            )
            return Nu_fd

        @staticmethod
        def average_nusselt_number(
            a: Quantity,
            b: Quantity,
            Pr: float,
            Gz: float
        ) -> float:
            """Average Nusselt number for simultaneously developing flow.
            No specific correlation available. Uses the correlation for the
            average Nusselt number of a circular tube as an approximation (but
            based on the hydraulic diameter of a rectangular tube).
            """
            UHF = RectangularTube.UniformHeatFlux
            CT_UHF = CircularTube.UniformHeatFlux
            Nu_fd = UHF.fully_developed_local_nusselt_number(a, b)
            Nu_fd_ct = CT_UHF.fully_developed_local_nusselt_number()
            Nu_avg = CT_UHF.average_nusselt_number(Pr, Gz)
            Nu_avg -= Nu_fd_ct  # subtract the circular tube's Nu_fd from its Nu_avg
            Nu_avg += Nu_fd  # and add the Nu_fd of the rectangular tube instead
            return Nu_avg

    class UniformWallTemperature:

        @staticmethod
        def fully_developed_local_nusselt_number(a: Quantity, b: Quantity) -> float:
            a = a.to('m').m
            b = b.to('m').m
            ar = min(a, b) / max(a, b)  # aspect ratio
            Nu_fd = 7.541 * (
                1
                - 2.610 * ar
                + 4.970 * (ar ** 2)
                - 5.119 * (ar ** 3)
                + 2.702 * (ar ** 4)
                - 0.548 * (ar ** 5)
            )
            return Nu_fd

        @staticmethod
        def average_nusselt_number(
            a: Quantity,
            b: Quantity,
            Pr: float,
            Gz: float
        ) -> float:
            """Average Nusselt number for simultaneously developing flow.
            No specific correlation available. Uses the correlation for the
            average Nusselt number of a circular tube as an approximation (but
            based on the hydraulic diameter of a rectangular tube).
            """
            UWT = RectangularTube.UniformWallTemperature
            CT_UWT = CircularTube.UniformWallTemperature
            Nu_fd = UWT.fully_developed_local_nusselt_number(a, b)
            Nu_fd_ct = CT_UWT.fully_developed_local_nusselt_number()
            Nu_avg = CT_UWT.average_nusselt_number(Pr, Gz)
            Nu_avg -= Nu_fd_ct  # subtract the circular tube's Nu_fd from its Nu_avg
            Nu_avg += Nu_fd  # and add the Nu_fd of the rectangular tube instead
            return Nu_avg


class AnnularDuct:

    class UniformHeatFlux:
        # No implementation available.
        pass

    class UniformWallTemperature:
        table = {
            'ri / ro': [0, 0.05, 0.10, 0.25, 0.50, 1.0],
            'Nu_fd_i': [float('nan'), 17.46, 11.56, 7.37, 5.74, 4.86],
            'Nu_fd_o': [3.66, 4.06, 4.11, 4.23, 4.43, 4.86]
        }
        interp_Nu_fd_i = interp1d(
            table['ri / ro'][1:], table['Nu_fd_i'][1:],
            kind='slinear'
        )
        interp_Nu_fd_o = interp1d(
            table['ri / ro'], table['Nu_fd_i'],
            kind='slinear'
        )

        @staticmethod
        def fully_developed_local_nusselt_number(ri: Quantity, ro: Quantity) -> float:
            """Nusselt number for fully developed laminar flow in circular tube
            annulus with outer surface insulated and inner surface at constant
            temperature.

            Uses the data from Incropera et al. FUNDAMENTALS OF HEAT AND MASS
            TRANSFER (7th Ed), p. 554, table 8.2 combined with interpolation.
            """
            ri = ri.to('m').m
            ro = ro.to('m').m
            ri_on_ro = ri / ro
            Nu_fd_i = AnnularDuct.UniformWallTemperature.interp_Nu_fd_i(ri_on_ro)
            return Nu_fd_i

        @staticmethod
        def average_nusselt_number(
            ri: Quantity,
            ro: Quantity,
            Pr: float,
            Gz: float
        ) -> float:
            """
            Average Nusselt number for laminar flow in circular tube
            annulus with outer surface insulated and inner surface at constant
            temperature.

            No specific correlation available. Uses the correlation for the
            average Nusselt number of a circular tube as an approximation (but
            based on the hydraulic diameter of an annular duct).
            """
            UWT = AnnularDuct.UniformWallTemperature
            CT_UWT = CircularTube.UniformWallTemperature
            Nu_fd_i = UWT.fully_developed_local_nusselt_number(ri, ro)
            Nu_fd_ct = CT_UWT.fully_developed_local_nusselt_number()
            Nu_avg = CT_UWT.average_nusselt_number(Pr, Gz)
            Nu_avg -= Nu_fd_ct  # subtract the circular tube's Nu_fd from its Nu_avg
            Nu_avg += Nu_fd_i  # and add the Nu_fd of the annular duct instead
            return Nu_avg
