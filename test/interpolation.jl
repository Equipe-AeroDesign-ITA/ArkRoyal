@assert begin
    ipts = [
        1.0 0.0;
        0.75 0.03;
        0.5 0.05;
        0.3 0.046;
        0.1 0.03;
        0.03 0.018;
        0.0 0.0;
        0.03 (-0.018);
        0.1 (-0.03);
        0.3 (-0.046);
        0.5 (-0.05);
        0.75 (-0.03);
        1.0 0.0;
    ]

    pts = Bezier_airfoil(
        ipts;
        rotate = false,
        normalize = false,
    )

    writedlm("interpolated.dat", pts)
    writedlm("interpolation_points.dat", ipts)

    true
end
