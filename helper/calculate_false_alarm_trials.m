function FA = calculate_false_alarm_trials(session_data)

FA = [];
CSminus_trials = find(session_data.pd.CSminus_allMask);
for t = CSminus_trials
    t_data = session_data.s(t);
    licks = sort([t_data.licksL t_data.licksR]) - t_data.CSon;
    if isempty(licks)
        FA = [FA 0];
    else
        licks = licks(1); % first lick counts
        if licks <= 1500 % 1500ms is the response window for rewards; false alarm
            FA = [FA 1];
        else
            FA = [FA 0];
        end
    end
end