function [tb, effectiveEmissivity, effectiveReflectivity] = ...
        radiative_stack(incidenceDeg, polarization, frequencyHz, ...
        layerThickness, layerTemperatureK, layerPermittivity, layerKind, ...
        bottomTemperatureK, bottomPermittivity, boundary, ...
        interfacePermittivityMode)
%RADIATIVE_STACK 保留上层下行源项的一次反射敏感性方案。

    sinTheta0 = sind(incidenceDeg);
    k0 = 2 * pi * frequencyHz / 299792458;
    layerCount = numel(layerThickness);

    effectiveSource = bottomTemperatureK;
    effectiveReflectivity = 0;
    lowerPermittivity = bottomPermittivity;
    lowerKind = 3;

    for layer = layerCount:-1:1
        currentPermittivity = layerPermittivity(layer);
        interfaceReflectivity = bottom_reflectivity( ...
            currentPermittivity, layerKind(layer), ...
            lowerPermittivity, lowerKind, boundary, ...
            interfacePermittivityMode, sinTheta0, polarization);

        effectiveSource = (1 - interfaceReflectivity) * ...
            effectiveSource;
        effectiveReflectivity = interfaceReflectivity + ...
            (1 - interfaceReflectivity)^2 * effectiveReflectivity;

        verticalIndex = sqrt(currentPermittivity - sinTheta0^2);
        transmission = exp(-2 * k0 * abs(imag(verticalIndex)) * ...
            layerThickness(layer));
        layerSource = layerTemperatureK(layer) * (1 - transmission);
        effectiveSource = layerSource * ...
            (1 + transmission * effectiveReflectivity) + ...
            transmission * effectiveSource;
        effectiveReflectivity = transmission^2 * effectiveReflectivity;

        lowerPermittivity = currentPermittivity;
        lowerKind = layerKind(layer);
    end

    topPermittivity = lowerPermittivity;
    if layerCount > 0 && string(interfacePermittivityMode) == "boundary"
        topPermittivity = boundary.top;
    end
    topReflectivity = fresnel_power( ...
        1, topPermittivity, sinTheta0, polarization);
    tb = (1 - topReflectivity) * effectiveSource;
    effectiveReflectivity = topReflectivity + ...
        (1 - topReflectivity)^2 * effectiveReflectivity;
    effectiveEmissivity = 1 - effectiveReflectivity;
    effectiveEmissivity = min(max(real(effectiveEmissivity), 0), 1);
end

function reflectivity = bottom_reflectivity( ...
        upperPermittivity, upperKind, lowerPermittivity, lowerKind, ...
        boundary, interfacePermittivityMode, sinTheta0, polarization)
    if upperKind == lowerKind
        reflectivity = 0;
    elseif string(interfacePermittivityMode) == "boundary" && ...
            upperKind == 2 && lowerKind == 3
        reflectivity = fresnel_power( ...
            boundary.iceBottom, lowerPermittivity, sinTheta0, polarization);
    elseif string(interfacePermittivityMode) == "boundary" && ...
            upperKind == 1 && lowerKind == 2
        reflectivity = fresnel_power( ...
            boundary.snowBottom, boundary.iceTop, sinTheta0, polarization);
    else
        reflectivity = fresnel_power( ...
            upperPermittivity, lowerPermittivity, sinTheta0, polarization);
    end
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
