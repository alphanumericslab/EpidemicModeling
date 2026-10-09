# Historical data

These files are copied unchanged from the attached repository's
`xprize-sample-data` folder. `OxCGRT_latest.csv` is a historical snapshot despite
its filename. The notebooks display its actual date range. Population counts
are from the bundled `Population2020` table.

- Oxford COVID-19 Government Response Tracker: https://github.com/OxCGRT/covid-policy-tracker
- XPRIZE Pandemic Response Challenge: https://www.xprize.org/challenge/pandemicresponse
- Johns Hopkins-format CSVs are supported by `read_covid19_data`; those upstream
  wide tables were not part of the supplied repository and are not required
  to run the examples.

Blank RegionName means the national entry. Preserve original NPI column names
when reading CSVs. Do not merge national and subnational rows before fitting.
Sample costs, forecasts, future intervention plans, and geographic lists are
research fixtures. Keep their upstream attribution and applicable terms.
