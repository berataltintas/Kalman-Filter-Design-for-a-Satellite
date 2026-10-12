function in_eclipse = check_eclipse(sat_pos, sun_vec)
  re = 6371e3;
  r_sun = dot(sat_pos, sun_vec);
  if r_sun > 0
    in_eclipse = 0;%aydınlıkta
  else
    r_perp = norm(sat_pos - r_sun * sun_vec);
    if r_perp < re
      in_eclipse = 1; %gölgede
    else
      in_eclipse = 0; %gece ama gölge konisi dışında
    end 
  end
end