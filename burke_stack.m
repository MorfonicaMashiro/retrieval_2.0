function [tb, effectiveEmissivity] = burke_stack( ...
        incidenceDeg, polarization, frequencyHz, layerThickness, ...
        layerTemperatureK, layerPermittivity, layerKind, bottomTemperatureK, ...
        bottomPermittivity, boundary, interfacePermittivityMode)
%BURKE_STACK Burke 非相干分层辐射的一阶递推。

    sinTheta0 = sind(incidenceDeg);
    k0 = 2 * pi * frequencyHz / 299792458;
    upward = bottomTemperatureK;
    upwardUnit = 1;
    lowerPermittivity = bottomPermittivity;
    lowerKind = 3;

    for layer = numel(layerThickness):-1:1
        currentPermittivity = layerPermittivity(layer);
        if layerKind(layer) == lowerKind
            reflectionBottom = 0;
        elseif string(interfacePermittivityMode) == "boundary" && ...
                layerKind(layer) == 2 && lowerKind == 3
            reflectionBottom = fresnel_power( ...
                boundary.iceBottom, bottomPermittivity, ...
                sinTheta0, polarization);
        elseif string(interfacePermittivityMode) == "boundary" && ...
                layerKind(layer) == 1 && lowerKind == 2
            reflectionBottom = fresnel_power( ...
                boundary.snowBottom, boundary.iceTop, ...
                sinTheta0, polarization);
        else
            reflectionBottom = fresnel_power( ...
                currentPermittivity, lowerPermittivity, ...
                sinTheta0, polarization);
        end

        verticalIndex = sqrt(currentPermittivity - sinTheta0^2);
        transmission = exp(-2 * k0 * abs(imag(verticalIndex)) * ...
            layerThickness(layer));
        sourceFactor = (1 - transmission) * ...
            (1 + reflectionBottom * transmission);
        upward = layerTemperatureK(layer) * sourceFactor + ...
            (1 - reflectionBottom) * upward * transmission;
        upwardUnit = sourceFactor + ...
            (1 - reflectionBottom) * upwardUnit * transmission;
        lowerPermittivity = currentPermittivity;
        lowerKind = layerKind(layer);
    end

    topPermittivity = lowerPermittivity;
    if ~isempty(layerThickness) && ...
            string(interfacePermittivityMode) == "boundary"
        topPermittivity = boundary.top;
    end
    reflectionTop = fresnel_power( ...
        1, topPermittivity, sinTheta0, polarization);
    tb = (1 - reflectionTop) * upward;
    effectiveEmissivity = (1 - reflectionTop) * upwardUnit;
    effectiveEmissivity = min(max(real(effectiveEmissivity), 0), 1);
end

function reflectivity = fresnel_power( ...
        upperPermittivity, lowerPermittivity, sinTheta0, polarization)
    qUpper = sqrt(upperPermittivity - sinTheta0^2);
    qLower = sqrt(lowerPermittivity - sinTheta0^2);
    if polarization == 1
        coefficient = (qUpper - qLower) / (qUpper + qLower);
    else
        coefficient = (lowerPermittivity * qUpper - ...
            upperPermittivity * qLower) / ...
            (lowerPermittivity * qUpper + upperPermittivity * qLower);
    end
    reflectivity = min(max(abs(coefficient)^2, 0), 1);
end
