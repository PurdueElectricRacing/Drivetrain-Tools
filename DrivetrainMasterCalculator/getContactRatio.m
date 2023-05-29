function[contact_ratio] = getContactRatio(module,pressureAngle,pinionTeeths,gearTeeths)
    addendum_radii_pinion = module + module.*pinionTeeths/2;
    base_radii_pinion = cosd(pressureAngle) * module.*pinionTeeths/2;
    lpinion = sqrt(addendum_radii_pinion.^2-base_radii_pinion.^2);
    
    addendum_radii_gear = module + module.* gearTeeths/2;
    base_radii_gear = cosd(pressureAngle) * module.*gearTeeths/2;
    lgear = sqrt(addendum_radii_gear.^2-base_radii_gear.^2);
    
    cdist = module.*(pinionTeeths+gearTeeths)/2;
    
    base_pitch = pi*module;
    contact_ratio = (lpinion+lgear-sind(pressureAngle)*cdist)./base_pitch;

end