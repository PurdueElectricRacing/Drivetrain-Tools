function[mp] = getContactRatio(module,pressureAngle,gear,pinion,cdist)
    dog = 2 * sqrt((cdist*sind(pressureAngle))^2+(gear.pitchDiameter/2*cosd(pressureAngle))^2);
    dop = 2 * sqrt((cdist*sind(pressureAngle))^2+(pinion.pitchDiameter/2*cosd(pressureAngle))^2);

    rog = dog/2;
    rop = dop/2;
    
    base_pitch = pi*module*cosd(pressureAngle);

    mp = (sqrt(rop^2-(pinion.pitchDiameter/2*cosd(pressureAngle))) + sqrt(rog^2-(gear.pitchDiameter/2*cosd(pressureAngle))) - cdist * sind(pressureAngle))/base_pitch;

    c5 = (rop^2-(pinion.pitchDiameter/2*cosd(pressureAngle))^2)^0.5;
    c6 = cdist * sind(pressureAngle);
    c1 = c6 - (rog^2 - (pinion.pitchDiameter/2*cosd(pressureAngle)))^0.5;
    z = c5- c1;
    mp = z/base_pitch;

end