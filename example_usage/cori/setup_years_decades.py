#!/usr/bin/env python
import os
import argparse
import numpy


def replace(inFileName, outFileName, replacements):
    with open(inFileName, 'rt') as fin:
        with open(outFileName, 'wt') as fout:
            for line in fin:
                for orig in replacements:
                    line = line.replace(orig, replacements[orig])
                fout.write(line)


def setup_years(start, end, model, scenario, ensemble):
    for tIndex, year in enumerate(range(start, end+1)):
        yearString = '{:04d}'.format(year)
        print(yearString)
        try:
            os.makedirs(yearString)
        except OSError:
            pass

        replacements = {'@tIndex': '{}'.format(tIndex),
                        '@modelLower': model.lower(),
                        '@modelUpper': model,
                        '@start': '{}'.format(start),
                        '@end': '{}'.format(end),
                        '@ensemble': ensemble,
                        '@scenario': scenario,
                        '@year': yearString}
        templateFileName = '../templates/config.template'
        outFileName = '{}/config.{}'.format(yearString, model.lower())
        replace(templateFileName, outFileName, replacements)
        templateFileName = '../templates/job_script_template.bash'
        outFileName = '{}/job_script.bash'.format(yearString)
        replace(templateFileName, outFileName, replacements)


def setup_decades(start, end, model, climStart=1995, climEnd=2014,
                  block=20):

    climOutFolder = '{:04d}-{:04d}'.format(climStart, climEnd)
    climFolders = ', '.join(['{:04d}'.format(year) for year in
                             range(climStart, climEnd+1)])

    histStart = numpy.array(range(start, climStart-1, block), dtype=int)
    histEnd = numpy.minimum(histStart+block-1, climStart-1)

    scenarioStart = numpy.array(range(climStart, end, block), dtype=int)
    scenarioEnd = numpy.minimum(scenarioStart+block-1, end)

    firstYears = numpy.append(histStart, scenarioStart)
    lastYears = numpy.append(histEnd, scenarioEnd)

    for firstYear, lastYear in zip(firstYears, lastYears):
        outFolder = '{:04d}-{:04d}'.format(firstYear, lastYear)
        folders = ', '.join(['{:04d}'.format(year) for year in
                             range(firstYear, lastYear+1)])
        print(outFolder)
        try:
            os.makedirs(outFolder)
        except OSError:
            pass

        replacements = {'@outFolder': outFolder,
                        '@folders': folders,
                        '@modelLower': model.lower(),
                        '@modelUpper': model,
                        '@years': outFolder,
                        '@tIndexMin': '{}'.format(firstYear-1850),
                        '@tIndexMax': '{}'.format(lastYear-1850),
                        '@climDecades': climOutFolder,
                        '@climFolders': climFolders,
                        '@climFirstTIndex': '0',
                        '@climLastTIndex': '{}'.format(climEnd-climStart)}
        templateFileName = '../templates/config.combine_template'
        outFileName = '{}/config.{}'.format(outFolder, model.lower())
        replace(templateFileName, outFileName, replacements)
        templateFileName = '../templates/job_script_combine_template.bash'
        outFileName = '{}/job_script.bash'.format(outFolder)
        replace(templateFileName, outFileName, replacements)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest="start", type=int,
                        help="start year of the time series")
    parser.add_argument("-e", dest="end", type=int,
                        help="end year of the time series")
    parser.add_argument("-m", dest="model", type=str,
                        help="name of the model")
    parser.add_argument("--scenario", dest="scenario", type=str,
                        help="scenario (rcp85, ssp585, etc.)")
    parser.add_argument("--ensemble", dest="ensemble", type=str,
                        help="ensemble member (r1i1p1, etc.)")
    args = parser.parse_args()

    setup_years(args.start, args.end, args.model, args.scenario, args.ensemble)
    setup_decades(args.start, args.end, args.model)
