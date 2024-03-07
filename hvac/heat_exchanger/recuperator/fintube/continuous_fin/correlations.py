import math


class PlainContinuousFinStaggeredTubeBank:
    
    @staticmethod
    def colburn_j_factor(
        Re: float,
        d_r: float,
        d_h: float,
        p_trv: float,
        p_lon: float,
        n_fin: float,
        n_rows: float
    ) -> float:
        """
        The Colburn j-factor correlation of Wang and Chi (2000) valid for plain 
        flat fins on staggered tube banks.

        Parameters
        ----------
        Re:
            Reynolds number based on the collar diameter of the fin.
        d_r:
            The collar diameter of the fin [m].
        d_h:
            The hydraulic diameter [m].
        p_trv:
            The lateral or transverse pitch, i.e., the spacing between tubes of 
            the same row [m].
        p_lon:
            The longitudinal pitch, i.e., the spacing between tubes of two 
            adjacent tube rows [m].
        n_fin:
            The number of fins per unit length = fin density [1/m]. It is the
            inverse of the fin spacing or pitch.
        n_rows:
            The number of rows.

        Notes
        -----
        Valid ranges of the parameters are:
        - Re: 300-20_000
        - collar diameter: 6.9-13.6 mm
        - hydraulic diameter: 1.30-9.37 mm
        - transverse pitch: 20.4-31.8 mm
        - longitudinal pitch: 12.7-32 mm
        - fin spacing: 1.0-8.7 mm
        - number of rows: 1-6

        Returns
        -------
        Colburn j-factor (float).
        """
        n_r = min(n_rows, 6)
        p_f = 1 / n_fin
        ln_Re = math.log(Re)
        if n_r < 2:
            c1 = 1.9 - 0.23 * ln_Re
            c2 = -0.236 + 0.126 * ln_Re
            k1 = (p_trv / p_lon) ** c1
            k2 = (p_f / d_r) ** -1.084
            k3 = (p_f / d_h) ** -0.786
            k4 = (p_f / p_trv) ** c2
            j = 0.108 * (Re ** -0.29) * k1 * k2 * k3 * k4
        else:
            c3 = (
                    -0.361 - (0.042 * n_r / ln_Re)
                    + 0.158 * math.log(n_r * (p_f / d_r) ** 0.41)
            )
            c4 = -1.224 - (0.076 * (p_lon / d_h) ** 1.42) / ln_Re
            c5 = -0.083 + 0.058 * n_r / ln_Re
            c6 = -5.735 + 1.21 * math.log(Re / n_r)
            k5 = n_r ** c4
            k6 = (p_f / d_r) ** c5
            k7 = (p_f / d_h) ** c6
            k8 = (p_f / p_trv) ** -0.93
            j = 0.086 * (Re ** c3) * k5 * k6 * k7 * k8
        return j

    @staticmethod
    def fanning_friction_factor(
        Re: float,
        d_r: float,
        p_trv: float,
        p_lon: float,
        n_fin: float,
        n_rows: float
    ) -> float:
        """Friction factor correlation by Wang and Chi (2000) for plain flat
        fins on staggered tube banks.

        Parameters
        ----------
        Re:
            Reynolds number based on collar diameter.
        d_r:
            Collar diameter [m].
        p_trv:
            Lateral or transverse pitch, i.e., spacing between tubes of the
            same row [m].
        p_lon:
            Longitudinal pitch, i.e., spacing between tubes of two adjacent tube
            rows [m].
        n_fin:
            The number of fins per unit length [1/m]
        n_rows:
            Number of rows.

        Notes
        -----
        Valid ranges of parameters:
        - Re: 300-20_000
        - collar diameter: 6.9-13.6 mm
        - hydraulic diameter: 1.30-9.37 mm
        - transverse pitch: 20.4-31.8 mm
        - longitudinal pitch: 12.7-32 mm
        - fin spacing: 1.0-8.7 mm
        - number of rows: 1-60

        Returns
        -------
        Fanning friction factor.
        """
        n_rows = min(n_rows, 6)
        p_f = 1 / n_fin  # fin pitch
        ln_Re = math.log(Re)
        c1 = -0.764 + 0.739 * (p_trv / p_lon) + 0.177 * (p_f / d_r) - 0.00758 / n_rows
        c2 = -15.689 + 64.021 / ln_Re
        c3 = 1.696 - 15.695 / ln_Re
        k1 = (p_trv / p_lon) ** c2
        k2 = (p_f / d_r) ** c3
        f = 0.0267 * (Re ** c1) * k1 * k2
        return f
