function ylabelunits(settings)

switch settings.units
    case 'prcchange'
        ylabel('% change from prestim')
    case 'log'
        ylabel('10*log10 unit change')
    case 'raw'
        ylabel('Change from prestim (respective units)')
    case 'zscore'
        ylabel('Normalized change from prestim')
end

end