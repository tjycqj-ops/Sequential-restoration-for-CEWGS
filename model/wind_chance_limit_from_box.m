function [wind_lower, wind_mean, wind_sigma, z_chance] = ...
    wind_chance_limit_from_box(point_forecast, error_lower, error_upper, ...
    epsilon_wind, box_coverage)
%WIND_CHANCE_LIMIT_FROM_BOX Deterministic lower limit for wind availability.
%   The supplied box is interpreted as a central prediction interval with
%   probability BOX_COVERAGE. A Gaussian forecast-error model is calibrated
%   from the interval midpoint and half-width. WIND_LOWER is the one-sided
%   EPSILON_WIND quantile used by the deterministic chance constraint.

    validateattributes(point_forecast, {'numeric'}, {'2d','real','finite','nonnegative'});
    validateattributes(error_lower, {'numeric'}, {'2d','real','finite'});
    validateattributes(error_upper, {'numeric'}, {'2d','real','finite'});
    validateattributes(epsilon_wind, {'numeric'}, {'scalar','real','finite','>',0,'<',0.5});
    validateattributes(box_coverage, {'numeric'}, {'scalar','real','finite','>',0.5,'<',1});

    if isrow(error_lower) && size(point_forecast,1) > 1
        error_lower = repmat(error_lower,size(point_forecast,1),1);
    end
    if isrow(error_upper) && size(point_forecast,1) > 1
        error_upper = repmat(error_upper,size(point_forecast,1),1);
    end
    if ~isequal(size(point_forecast),size(error_lower),size(error_upper))
        error('WindChance:DimensionMismatch', ...
            'Forecast and expanded error-bound arrays must have identical sizes.');
    end
    if any(error_upper(:) <= error_lower(:))
        error('WindChance:InvalidBox', ...
            'Every upper forecast-error bound must exceed its lower bound.');
    end

    error_center = (error_upper+error_lower)/2;
    error_halfwidth = (error_upper-error_lower)/2;
    box_tail_probability = (1-box_coverage)/2;
    z_box = -sqrt(2)*erfcinv(2*(1-box_tail_probability));
    z_chance = -sqrt(2)*erfcinv(2*(1-epsilon_wind));

    wind_sigma = error_halfwidth/z_box;
    wind_mean = point_forecast+error_center;
    wind_lower = wind_mean-z_chance*wind_sigma;
    wind_lower = min(max(wind_lower,0),1);
end


