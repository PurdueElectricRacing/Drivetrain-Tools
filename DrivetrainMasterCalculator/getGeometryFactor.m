function[J] = getGeometryFactor(N)
    J = -0.1112*(log(N)).^2+0.5533*(log(N))-0.2126;
end