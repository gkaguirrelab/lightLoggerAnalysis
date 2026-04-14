function result = msTimestampSlidingWindow(ms_t, ms_v, window_size_ns, operation)
    arguments
        ms_t (:,1) double
        ms_v (:,1) double
        window_size_ns (1,1) double {mustBePositive}
        operation (1,1) string {mustBeMember(operation, ["mean", "median", "std"])}
    end

    % Define the operator
    % we will use for the window 
    switch operation
        case "mean"
            operator = @(x) mean(x, "omitnan");
        case "median"
            operator = @(x) median(x, "omitnan");
        case "std"
            operator = @(x) std(x, "omitnan");
        otherwise
            error("Unsupported operation: %s", operation);
    end

    % Initialize output
    n_readings = numel(ms_t);
    result = nan(n_readings, 1);

    % Any given value i should be in the center 
    % of the window, so this tells us 
    % the approximate left and right side of the window
    half_window_ns = window_size_ns / 2;

    % Bounds pointers 
    l = 1;
    r = 1;

    % Move across the values in the readings
    for ii = 1:n_readings
        center_t = ms_t(ii);
        left_limit = center_t - half_window_ns;
        right_limit = center_t + half_window_ns;

        % Move left bound rightward until it is inside the window
        while l <= n_readings && ms_t(l) < left_limit
            l = l + 1;
        end

        % Make sure right bound is at least i
        if r < i
            r = i;
        end

        % Move right bound rightward while still inside the window
        while r < n_readings && ms_t(r + 1) <= right_limit
            r = r + 1;
        end

        % Apply operation over centered timestamp window
        result(ii) = operator(ms_v(l:r));
    end
end