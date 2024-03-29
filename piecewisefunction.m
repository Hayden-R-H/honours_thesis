function y = piecewisefunction(t)
    t1 = (t >= 0 & t < 1 / 4);
    t2 = (t >= 1 / 4 & t < 1 / 2);
    t3 = (t >= 1 / 2 & t <= 3 / 4);
    t4 = (t >= 3 / 4 & t < 1);
    y(t1) = 2 * t(t1);
    y(t2) = (31 / 8) * (t(t2) - 1/4);
    y(t3) = (31 / 8) * (t(t3) - 3 / 4) + 1;
    y(t4) = 2 * (t(t4) - 1) + 1;