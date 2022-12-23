function [min_sunteeth,rim_thickness,effectiveInternalRadius] = minSunTeeth(pressureAngle,module,sun,motor)
    backup_ratio = 1.2;
    deddendum_factor = sun.whole_depth_factor - sun.addendum_factor;
    tooth_height = sun.whole_depth_factor * module;
    rim_thickness = tooth_height * backup_ratio;
    effectiveInternalRadius = sqrt((sun.internalDiameter / 2 + motor.keyProjHeight)^2 + (motor.keyWidth / 2)^2);
    pitchDiameter = 2 * (rim_thickness + deddendum_factor * module + 1000 * effectiveInternalRadius);
    min_sunteeth_rim = ceil(pitchDiameter / module);
    min_sunteeth_undercut = ceil(2 / (sind(pressureAngle))^2);
    min_sunteeth = max(min_sunteeth_rim,min_sunteeth_undercut);
end