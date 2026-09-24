# DigiSizer Parser

`digisizer_parser.py` converts a Micromeritics Saturn DigiSizer legacy Excel export (`.XLS`) into a tidy CSV file for analysis.

The script extracts the sample identifier and the particle-size distribution table from the instrument export. It writes a CSV containing the particle diameter, volume-frequency distribution and cumulative finer volume percent.

## Files supplied

| File | Purpose |
| --- | --- |
| `T1SG1 - 005-941.pdf` | Readable instrument report for sample T1SG1, including analysis settings, distribution summaries and plots. |
| `005-941(1).XLS` | Raw legacy Excel export containing the particle-size data. |
| `digisizer_parser(1).py` | Python script that extracts the data from the XLS export into CSV format. |

## Requirements

- Python 3.8 or later
- [`pandas`](https://pandas.pydata.org/)
- [`xlrd`](https://pypi.org/project/xlrd/), which reads the older `.XLS` format used by the instrument export

Check that Python is available:

```bash
python --version
```

Install the dependencies:

```bash
python -m pip install pandas xlrd
```

On some Linux and macOS systems, use `python3` instead of `python`. If you use a virtual environment, activate it before installing packages and running the script.

## Run the parser

Place the script and the `.XLS` file in the same folder, then open a terminal in that folder.

### Default output name

```bash
python "digisizer_parser(1).py" "005-941(1).XLS"
```

The script writes `005-941(1).csv` in the same folder as the source XLS file.

### Choose an output name

```bash
python "digisizer_parser(1).py" "005-941(1).XLS" --output T1SG1_particle_size.csv
```

The short option is equivalent:

```bash
python "digisizer_parser(1).py" "005-941(1).XLS" -o T1SG1_particle_size.csv
```

Quote filenames that contain spaces or parentheses, as in the examples above.

## Output

The output CSV includes:

- `sample_id`: sample name extracted from the XLS file
- particle diameter and volume frequency percent
- particle diameter and cumulative finer volume percent

For the supplied file, `sample_id` is `T1SG1`. The cumulative finer volume percent should increase from approximately 0 to 100.

## How it works

The script:

1. Reads the XLS file without assuming a fixed header row.
2. Finds the `Sample:` field to obtain the sample identifier.
3. Finds the row beginning with `Particle Diameter (µm)`.
4. Extracts the distribution table until the first blank row.
5. Adds `sample_id` and writes the results as CSV.

## Troubleshooting

| Problem | Suggested action |
| --- | --- |
| `python` is not recognised | Install Python, or try `python3` instead. |
| `ModuleNotFoundError` for `pandas` or `xlrd` | Run `python -m pip install pandas xlrd`, then repeat the command. |
| File not found | Check the filename, current folder and quotation marks. Use a full file path if necessary. |
| The script cannot find `Sample:` or `Particle Diameter (µm)` | Check that the XLS export follows the standard DigiSizer layout expected by the script. A differently structured export may need the parser to be adapted. |

## Use with other DigiSizer exports

Use the same command pattern with any compatible DigiSizer `.XLS` export. The script derives `sample_id` from the file contents, so no manual sample-name entry is required. Keep the original XLS and PDF report as source records; use the CSV as the analysis-ready file.
