export BSpline_airfoil

"""
```
    function BSpline_airfoil(
        pts::AbstractMatrix;
        discretization::Int = 250,
        rotate::Bool = true,
        normalize::Bool = true,
        close_TE::Bool = true,
    )
```

Interpolate airfoil from points in `(npts, 2)` matrix `pts` using cubic
B-Spline. 

Interpolates with equal number of grid points between interpolation points
(i. e. more spline control points concentrated in a given region imply a finer
local discretization).

Normalizes the airfoil to chord 1 (highly recommended) if `normalize` is set to `true`.
Aligns the camberline with the `x` axis (sets `α = 0`) if `rotate` is set to `true`.

Returns matrix of points to be fed to `Airfoil` constructor.
"""
function BSpline_airfoil(
    pts::AbstractMatrix;
    discretization::Int = 250,
    rotate::Bool = true,
    normalize::Bool = true,
    close_TE::Bool = true,
)

    xint = extrapolate(
        interpolate(pts[:, 1], BSpline(Cubic(Line(OnGrid())))),
        Flat()
    )
    yint = extrapolate(
        interpolate(pts[:, 2], BSpline(Cubic(Line(OnGrid())))),
        Line()
    )

    discretization = collect(
        LinRange(
            1, size(pts, 1), discretization
        )
    )

    npts = [
        xint.(discretization) yint.(discretization)
    ]

    if rotate
        iTE = argmax(npts[:, 1])
        iLE = argmin(npts[:, 1])

        θ = atan(
            npts[iLE, 2] - npts[iTE, 2], - npts[iLE, 1] + npts[iTE, 1]
        )

        R = [
            cos(θ) (-sin(θ));
            sin(θ) cos(θ)
        ]

        npts = (npts * R')
    end

    if normalize
        npts[1, :] .-= minimum(npts[1, :])

        npts ./= (maximum(npts[:, 1]) - minimum(npts[:, 1]))
    end

    if close_TE
        let m = (npts[1, :] .+ npts[end, :]) ./ 2
            npts[1, :] .= m
            npts[end, :] .= m
        end
    end

    npts

end
