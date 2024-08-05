function [min_sunteeth,rim_thickness,effectiveInternalRadius] = minSunTeeth(pressureAngle,module,sun,motor)
    backup_ratio = 1.2; % This is the generally accepted value for ratio of rim thickness to tooth height (AGMA standards)
    deddendum_factor = sun.whole_depth_factor - sun.addendum_factor; 
    tooth_height = sun.whole_depth_factor * module;
    rim_thickness = tooth_height * backup_ratio;

    % IMPORTANT - THIS CALCULATION ASSUMES THE KEYED SHAFT OF THE NOVA 15
    % IS PRESENT. UPDATE effectiveInternalRadius TO THE MAXIMUM EXTERNAL
    % DIAMETER OF THE SHAFT OF THE MOTOR BEING CONSIDERED.
    effectiveInternalRadius = sqrt((sun.internalDiameter / 2 + motor.keyProjHeight)^2 + (motor.keyWidth / 2)^2);

    % Calculate pitch diameter from the known shaft diameter, rim thickness
    % and tooth geometry.
    pitchDiameter = 2 * (rim_thickness + deddendum_factor * module + 1000 * effectiveInternalRadius);

    % Calculate the minimum number of teeth based on rim thickness
    min_sunteeth_rim = ceil(pitchDiameter / module);

    % Calculate the minimum number of teeth based on undercutting (lookup
    % formula online)
    min_sunteeth_undercut = ceil(2 / (sind(pressureAngle))^2);

    % Return the true minimum teeth as the maximum of the two calculated
    % values, as this would meet both requirements
    min_sunteeth = max(min_sunteeth_rim,min_sunteeth_undercut);
end