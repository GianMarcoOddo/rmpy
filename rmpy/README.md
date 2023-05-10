# Risk Management Python (rmpy) Package

The `rmpy` package is a comprehensive and powerful tool designed for risk management and quantitative finance in Python. It provides a suite of functionalities to perform essential risk assessments, calculations, and analyses on financial assets and portfolios. This package streamlines the process of evaluating the risk and performance of various assets, allowing users to make informed decisions for managing their investments.

---

## Installation

```python
pip install rmpy
```
---
## If the installation process encounters an error, make sure to have xcode installed by typing the following command in the terminal

```python
xcode-select --install
```
## Github Repository

## https://github.com/GianMarcoOddo/rmpy

---

## Key Features

- Parametric and Non-Parametric Value at Risk (pVaR) calculations for single assets and portfolios.
- Historical pVaR and npVaR calculation using Yahoo Finance data for single assets and portfolios.
- Confidence level pVaR and npVaR calculations for single assets and portfolios.
- Support for daily and weekly frequency data and customizable intervals.
- Calculation of marginal, component, and relative component Value at Risk (VaR) for portfolios.
- A user-friendly interface with the option to display or suppress output as needed.

## Who is it for?

The `rmpy` package is ideal for quantitative analysts, portfolio managers, risk managers, and anyone interested in financial risk management. The package is easy to use, well-documented, and compatible with various data sources and formats. It is a powerful addition to the Python ecosystem for risk management and a valuable resource for those looking to enhance their understanding of financial risk.

Get started with `rmpy` today and take control of your financial risk management!

# 1. NpVaR Module

This module calculates the Non-parametric Value-at-Risk (NpVaR) and associated functions for single assets or portfolios using historical returns. 
Functions with the `yf_` prefix use data from the Yahoo Finance API, whereas the others use the provided returns. The `_single` functions are for individual assets, and the `_port` functions are for portfolios. The `_conflevel_` functions calculate the NpVaR with a specified confidence level. The `_summary_` functions provide summaries for NpVaR calculations. The `_marg_NpVaRs_scale_factor` functions calculate the marginal NpVaRs given a scaling factor.
Please refer to the code below for the syntax and examples on how to use each function. The input parameters and their usage are described within the comments.

It contains the following functions:

`1. yf_npVaR_single`

- yf_npVaR_single(ticker, position, start_date, end_date, freq="daily", alpha=0.01, kind="abs", display=True)

This function yf_npVaR_single calculates the quantile non-parametric value at risk (VaR) for a SINGLE asset position using historical prices obtained from Yahoo Finance (yfinance).

Non-parametric VaR is a method of calculating the minimum amount of loss that an asset is likely to experience at a given confidence level, using historical data and without assuming any particular probability distribution.

#### Args:

- ticker: the asset symbol or ticker of the company for the npVaR will be calculated
- position: the size of the position (value).
- start_date: the starting date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object.
- end_date: the ending date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object. By default, this is set to "today" 
which will use the current date as the ending date.
- freq: the frequency of the historical prices. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".
- alpha: the significance level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
- kind: the type of VaR calculation to use. This can be set to "abs" for absolute VaR or "rel" for relative VaR. By default, this is set to "abs".
- display: a boolean value or string representing whether to display the calculated VaR. This can be set to True or False. By default, this is set to True.

### Example:

```python
from rmpy.NpVaR import yf_npVaR_single
ticker = ['AAPL'] 
position = [-1000] 
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
alpha = 0.01
kind = "abs"
display = True
yf_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha,kind = kind, display=display)
# OR
VaR = yf_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha,kind = kind, display=False)
print(VaR)
```
This example calculates the daily non-parametric VaR for a short position of 1000 in Apple Inc, with a confidence level of 99%.

`2. yf_npVaR_port`

- yf_npVaR_port(tickers, positions, start_date, end_date, freq="daily", alpha=0.01, kind="abs", display=display)

This function calculates the quantile non-parametric Value at Risk (VaR) for a portfolio of assets using historical prices obtained from Yahoo Finance (yfinance).

Non-parametric VaR is a way of estimating the potential downside risk of a portfolio that does not rely on any specific assumptions about the shape or parameters of the underlying probability distribution.

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portfolio. #### Note that all the TICKERS provided should be part of the portfolio whose VaR is being calculated ####
- positions: A list of integers or floats representing the positions of the assets in the portfolio. The length of this list should be the same as the 'tickers' list.
- start_date: A string representing the start date for the historical price data in the format 'YYYY-MM-DD'.
- end_date: A string representing the end date for the historical price data in the format 'YYYY-MM-DD'. By default, it is set to "today".
- freq: A string representing the frequency of the price data. By default, it is set to "daily".
- alpha: A float representing the confidence level for the VaR calculation. By default, it is set to 0.01.
- kind: A string representing the type of VaR calculation to perform. It can be either "abs" for absolute VaR or "rel" for relative VaR. By default, it is set to "abs".
- display: A boolean indicating whether to print the VaR results or not. By default, it is set to True.

### Example:

```python
from rmpy.NpVaR import yf_npVaR_port
tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
alpha = 0.01
kind = "abs"
display = True
yf_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, kind=kind, display=display)
# OR
VaR = yf_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, kind=kind, display=False)
print(VaR)
```
This example calculates the daily non-parametric VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

`3. npVaR_single`

- VaR = npVaR_single(returns, position, alpha=0.01, kind = "abs")

This function npVaR_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE asset using historical returns data (you can use every type of assets e.g stock, options, bonds, etc.).

Non-parametric VaR is a way to estimate the potential loss that an asset may experience over a given period, based on a ranking of the historical returns without making any assumptions about the underlying probability distribution.

#### Args:

- returns: a pandas Series or NumPy array containing the historical returns of the asset.
- position: the size of the position in units of the asset. 
- alpha: the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
- kind: the type of VaR calculation to use. This can be set to "abs" for absolute VaR or "rel" for relative VaR. By default, this is set to "abs".

### Example:

```python
from rmpy.NpVaR import npVaR_single
import numpy as np
returns = np.array([-0.02, 0.03, -0.01, 0.015, -0.002, 0.001, 0.008, 0.002, -0.006, 0.009]) # Replace this with actual returns
position = -1000
alpha = 0.01
kind = "abs"
VaR = npVaR_single(returns, position, alpha=alpha, kind=kind)
print(VaR)
```
This example calculates the absolute non-parametric VaR consisting of a single short position of 1000 for the given returns.

`4. npVaR_port`

- VaR = npVaR_port(returns, position, alpha=0.01, kind = "abs")

This function npVaR_port calculates the quantile non-parametric value at risk (VaR) for a PORTFOLIO of assets using historical returns data (you can use every type of assets e.g stock, options, bonds, etc.).

Non-parametric VaR is a risk management tool that uses a variety of techniques, such as bootstrapping, to estimate the distribution of returns based on historical data, without making any assumptions about the underlying probability distribution.

#### Args:

- returns: a pandas Series or NumPy array containing the historical returns of the portfolio.  #### note that all the RETURNS provided should be part of the portfolio whose VaR is being calculated ####
- positions: the size of the positionS in units of the portfolio. This can be a single value or an array of values corresponding to each element in the returns argument.
- alpha: the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
- kind: the type of VaR calculation to use. This can be set to "abs" for absolute VaR or "rel" for relative VaR. By default, this is set to "abs".

### Example:

```python
from rmpy.NpVaR import npVaR_port
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 2500, 7000]
alpha = 0.01
kind = "abs"
VaR = npVaR_port(returns, positions, alpha=alpha, kind=kind)
print(VaR)
```
This example calculates the relative non-parametric VaR consisting of a portfolio with short positions of 1000 in the first asset, long positions of 2500 in the second asset, and long positions of 7000 in the third asset.

`5. yf_conflevel_npVaR_single`

- yf_conflevel_npVaR_single(ticker, position, start_date, end_date, freq="daily", alpha= 0.01,confidence_level = 0.95, display=True)

This function yf_conflevel_npVaR_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE asset using Yahoo Finance data. It first downloads historical price data from Yahoo Finance, calculates the returns of the asset, and then calculates the VaR at its confidence level (lower and upper bound).

Non-parametric VaR is a risk management technique used to estimate the potential losses in an asset, based on empirical data rather than any assumptions about the underlying distribution of the returns.

#### Args:

- ticker: the symbol of the asset to calculate VaR for.
- position: the size of the position in units of the asset.
- start_date: the start date of the historical data to download. This should be a string in the format "YYYY-MM-DD".
- end_date: the end date of the historical data to download. This should be a string in the format "YYYY-MM-DD". By default, this is set to "today".
- freq: the frequency of the data to download. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".
- alpha: the level of significance for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
- confidence_level: the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability that the true VaR is within the calculated VaR interval. By default, this is set to 0.95 (95%).
- display: a boolean value indicating whether to display the VaR calculation result. By default, this is set to True.

### Example:

```python
from rmpy.NpVaR import yf_conflevel_npVaR_single
ticker = ['AAPL']
position = [2500]
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
alpha = 0.01
confidence_level = 0.95
display = True
yf_conflevel_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha, confidence_level=confidence_level, display=display)
# OR
VaR = yf_conflevel_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha, confidence_level=confidence_level, display=False)
print(VaR)
```
This example calculates the absolute non-parametric VaR with its lower and upper bound (95% confidence) for a position of 2500 in "AAPL".

`6. yf_conflevel_npVaR_port`

-  yf_conflevel_npVaR_port(tickers, positions, start_date, end_date, freq= "daily",alpha=0.01,confidence_level=0.95, display=display)

This function yf_conflevel_npVaR_port calculates the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using Yahoo Finance data. It first downloads historical price data from Yahoo Finance for each asset in the portfolio, calculates the returns of each asset, and then  calculates the portfolio VaR, the lower and the upper bound at a specified confidence level and alpha level.

Non-parametric VaR is a risk management approach that uses non-parametric methods such as kernel density estimation to estimate the distribution of returns without assuming any specific distribution.

#### Args:

- tickers: a list of symbols of the assets in the portfolio to calculate VaR for.
- positions: a list or array containing the sizes of the positions for each asset in the portfolio. These sizes are in units of the assets and can be positive or negative.
- start_date: the start date of the historical data to download. This should be a string in the format "YYYY-MM-DD".
- end_date: the end date of the historical data to download. This should be a string in the format "YYYY-MM-DD". By default, this is set to "today".
- freq: the frequency of the data to download. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".
- alpha: the level of significance for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
- onfidence_level: the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability that the true VaR is within the calculated VaR interval. By default, this is set to 0.95.
- display : a boolean value indicating whether to display the VaR calculation result. By default, this is set to True.

### Example:

```python
from rmpy.NpVaR import yf_conflevel_npVaR_port
tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
alpha = 0.01
confidence_level = 0.95
display = True
yf_conflevel_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, confidence_level=confidence_level, display=display)
# OR
VaR = yf_conflevel_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, confidence_level=confidence_level, display=False)
print(VaR)
```
This example calculates the daily non-parametric VaR with its lower and upper bound (95% confidence) for a Portfolio with short position of 1000 in Apple Inc and 5000 in Microsoft Corp.

`7. conflevel_npVaR_single`

- VaR = conflevel_npVaR_single(returns, position,confidence_level= 0.95, alpha=0.01)

This function conflevel_npVaR_single calculates the quantile non-parametric value at risk (VaR) for a SINGLE asset using historical returns data (you can use every type of assets e.g stock,options,bonds, ecc.) with a specific confidence interval and alpha value. It also caluculates the lower and upper bound for the non-parametric VaR.

Non-parametric VaR is a measure of the worst-case loss that an asset is likely to experience at a given confidence level, using empirical data to estimate the distribution of returns without relying on any particular parametric model.

#### Args:

 - returns: A NumPy array or pd.Series of historical returns of a single asset.
- position: The size of the position in the asset. If position is greater than 0, the function assumes the position is long, and if position is less than or equal to 0, the function assumes the position is short.
- confidence_level: The confidence level for the VaR calculation (default value is 0.95).
- alpha: The level of significance for the VaR calculation (default value is 0.01).

### Example:

```python
from rmpy.NpVaR import conflevel_npVaR_single
import numpy as np
returns = np.random.uniform(-0.03, 0.03, size=1000) # Replace this with actual returns
position = -1000
confidence_level = 0.95
alpha = 0.01
VaR = conflevel_npVaR_single(returns, position, confidence_level=confidence_level, alpha=alpha)
print(VaR)
# OR
VaR = yf_conflevel_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, confidence_level=confidence_level, display=False)
print(VaR)
```
This example calculates the non-parametric VaR and its lower and upper bound (95% confidence interval) consisting of a position of -1000 (short position) in the given asset.

`8. conflevel_npVaR_port`

- VaR = conflevel_npVaR_port(returns, positions ,confidence_level= 0.95, alpha=0.01)

This function conflevel_npVaR_port calculates the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using historical returns data (you can use every type of assets e.g stock, options, bonds, ecc.) with a specific confidence interval and alpha value. It also calculates the lower and the upper bound of this estimation. 

Non-parametric VaR is a way of estimating the potential downside risk of a portfolio that is based on quantiles of the historical returns rather than any specific parametric model.

#### Args:

- returns: A NumPy array of historical returns for the portfolio.
- positions: A NumPy array of positions for each asset in the portfolio.
- confidence_level: The confidence level for the VaR calculation (default value is 0.95).
- alpha: The level of significance for the VaR calculation (default value is 0.01).

### Example:

```python
from rmpy.NpVaR import conflevel_npVaR_port
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 2000, -3000]
confidence_level = 0.95
alpha = 0.01
VaR = conflevel_npVaR_port(returns, positions, confidence_level=confidence_level, alpha=alpha)
print(VaR)
```
This example calculates the non-parametric VaR consisting and its lower and upper bound (95% CI) of portfolio consisting of a postion of -1000 (short position) in the first asset, 2000 (long position) in the second asset and -3000 (short position) in the third one.

`9. yf_npVaR_summary_single`

- yf_npVaR_summary_single(ticker, position, start_date, end_date, freq="daily", alpha=0.01, display=display)
    
This function yf_npVaR_summary_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE stock position using historical prices obtained from Yahoo Finance (yfinance). This function is useful to visualize some risk metrics of the single asset npVaR.

Non-parametric VaR is a risk management technique that uses a variety of statistical techniques, such as historical simulation, to estimate the distribution of returns without relying on any specific probability distribution.

#### Args:

- ticker: The stock symbol or ticker of the company for which we want to calculate the VaR. This can be a string (for a single ticker) or a list or array of tickers (for multiple tickers).
- position: The size of the position in shares or currency. This can be a single value or an array of values corresponding to each ticker in the ticker argument.
- start_date: The starting date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object.
- end_date: The ending date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object. By default, this is set to "today" which will use the current date as the ending date.
- freq: The frequency of the historical prices. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".
- alpha: The confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
- display: A boolean value or string representing whether to display the calculated VaR. This can be set to True, "True", "T", False, "False", or "F". By default, this is set to True.

### Example:

```python
from rmpy.NpVaR import yf_npVaR_summary_single
ticker = 'AAPL'
position = -1000
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
alpha = 0.01
display = True
yf_npVaR_summary_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha, display=display)
# OR
VaR = yf_npVaR_summary_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha, display=False)
```
This example calculates several metrics for the non-parametric VaR summary for a short position of 1000 shares in Apple Inc with a confidence level of 99% using historical daily prices from January 1, 2020, to December 31, 2021, obtained from Yahoo Finance.

`10. yf_npVaR_summary_port`

- yf_npVaR_summary_port(tickers, positions, start_date, end_date, freq="daily", alpha=0.01, display=True)

 This function 'yf_npVaR_summary_port' is designed to calculate the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using historical prices obtained from Yahoo Finance (yfinance). This function is useful to visualize several risk metrics of the Portfolio npVaR.
 
Non-parametric VaR is a risk management technique that estimates the minimum loss a portfolio is likely to experience based on historical data, without assuming any specific probability distribution or parametric model.

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portfolio. #### note that all the TICKERS provided should be part of the portfolio whose VaR is being calculated ####
- positions: A list of integers or floats representing the positions of the assets in the portfolio. The length of this list should be the same as the 'tickers' list.
- start_date: A string representing the start date for the historical price data in the format 'YYYY-MM-DD'.
- end_date: A string representing the end date for the historical price data in the format 'YYYY-MM-DD'. By default, it is set to "today".
- freq: A string representing the frequency of the price data. By default, it is set to "daily".
- alpha: A float representing the confidence level for the VaR calculation. By default, it is set to 0.01.
- kind: A string representing the type of VaR calculation to perform. It can be either "abs" for absolute VaR or "rel" for relative VaR. By default, it is set to "abs".
- display: A boolean indicating whether to print the VaR results or not. By default, it is set to True.

### Example:

```python
from rmpy.NpVaR import yf_npVaR_summary_port
tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
alpha = 0.01
display = True
yf_npVaR_summary_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, display=display)
# OR
VaR = yf_npVaR_summary_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, display=False)
```
This example calculates several metrics for the non-parametric VaR for a portoflio consisting of a short position of 1000 in Apple Inc and a long position of 5000 in Miscrosoft using historical daily prices from January 1, 2020, to December 31, 2021, obtained from Yahoo Finance.

`11. npVaR_summary_single`

- summary = npVaR_summary_single(returns, position, alpha = 0.01)

This function called npVaR_summary_single that calculates the quantile Non-Parametric Value at Risk (VaR) and several risk metrics for a SINGLE asset and returns a summary of the VaR calculation.

Non-parametric VaR is a risk management technique that uses order statistics, such as percentiles, to estimate the potential losses that an asset may experience without making any assumptions about the underlying distribution.

#### Args:

- returns: A numpy array or pandas series of historical returns for the asset.
- position: The size of the position in the asset.
- alpha: The level of significance for the VaR calculation. Default value is 0.01.

### Example:

```python
import numpy as np
from rmpy.NpVaR import npVaR_summary_single
returns = np.array([-0.02, 0.03, -0.01, 0.015, -0.002, 0.001, 0.008, 0.002, -0.006, 0.009])
position = 10000
alpha = 0.01
summary = npVaR_summary_single(returns, position, alpha=alpha)
print(summary)
```
This example calculates all the metrics for the non-parametric VaR for a single asset with a given set of historical returns and position.

`12. npVaR_summary_port`

- summary = npVaR_summary_port(returns, positions, alpha=alpha)

The function npVaR_summary_port calculates the quantile non-parametric Value at Risk (VaR) and several risk metrics for a npVaR of a portfolio.

Non-parametric VaR is a risk management tool that uses a variety of techniques, to estimate the distribution of returns based on historical data, without making any assumptions about the underlying probability distribution.

#### Args:

- returns: A numpy array or pandas series of historical returns for the portfolio.
- position: which is a pandas DataFrame or numpy array containing the value held for each asset in the portfolio
- alpha: which is the confidence level for the VaR calculation (default value is 0.01).

### Example:

```python
import numpy as np
from rmpy import npVaR_summary_port
returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
positions = [1000, 500, -1500]
alpha = 0.01
summary = npVaR_summary_port(returns, positions, alpha=alpha)
print(summary)
```
This example shows how to use the npVaR_summary_port function with a 10x3 numpy array of historical returns and a list of positions for 3 assets in the portfolio. The function calculates the non-parametric Value at Risk (VaR) and other risk measures, returning them in a dictionary.

`13. yf_marg_NpVaRs_scale_factor`

- yf_marg_NpVaRs_scale_factor(tickers, positions, start_date, end_date, scale_factor=0.1, freq="daily", alpha=0.01, display=True)

The function yf_marg_NpVaRs_scale_factor calculates the Marginal Non-Parametric VaRs and their changes for a PORTFOLIO of assets over a given period of time, using a specified scale factor,using historical prices obtained from Yahoo Finance (yfinance). New position = old one + (old one * scale_factor) in the same direction (- if short and + if long).

Marginal Non-Parametric VaRs refer to the individual VaRs of each asset in a portfolio, calculated without assuming any specific probability distribution. They help identify which assets contribute the most to the overall risk of the portfolio and can be used to make informed decisions about risk management strategies, such as diversification or hedging.

#### Args:

- tickers: A list of tickers for the assets in the portfolio.
- positions: A list of positions for each asset in the portfolio.
- start_date: The start date of the analysis period.
- end_date: The end date of the analysis period (optional, default is "today").
- scale_factor: The scale factor used to calculate the Marginal Non-Parametric VaRs (optional, default is 0.1).
- freq: The frequency of the data (optional, default is "daily").
- alpha: The confidence level used in calculating the VaRs (optional, default is 0.01).
- display: A boolean indicating whether to print the resulting VaR changes (optional, default is True).

### Example:

```python
from rmpy.NpVaR import yf_marg_NpVaRs_scale_factor
tickers = ["AAPL", "MSFT", "GOOG"]
positions = [1000, 1500, 2000]
start_date = "2021-01-01"
end_date = "2021-12-31"
scale_factor = 0.1
freq = "daily"
alpha = 0.01
display = True
yf_marg_NpVaRs_scale_factor(tickers, positions, start_date, end_date, scale_factor=scale_factor, freq=freq, alpha=alpha, display=display)
# OR
NpVaR_changes = yf_marg_NpVaRs_scale_factor(tickers, positions, start_date, end_date, scale_factor=scale_factor, freq=freq, alpha=alpha, display=False)
print(NpVaR_changes)
```
This example demonstrates how to use the yf_marg_NpVaRs_scale_factor function with a list of three tickers, a list of positions for each asset, a start date, and an end date for the analysis period. The function calculates the Marginal Non-Parametric VaRs and their changes for the given portfolio using historical prices from Yahoo Finance.

`14. marg_NpVaRs_scale_factor`

- NpVaR_changes = marg_NpVaRs_scale_factor(returns, positions, alpha=0.01)

The function marg_NpVaRs_scale_factor calculates the marginal changes in the non-parametric value at risk (NpVaR) of a PORTFOLIO based on a specified scale factor. New position = old one + (old one * scale_factor) in the same direction (- if short and + if long).

Marginal Non-Parametric VaRs refer to the individual VaRs of each asset in a portfolio, calculated without assuming any specific probability distribution. They help identify which assets contribute the most to the overall risk of the portfolio and can be used to make informed decisions about risk management strategies, such as diversification or hedging.

#### Args:

- returns: A numpy array or pandas series of historical returns for the portfolio.
- positions: A list or array containing the current positions of the assets in the portfolio.
- scale_factor: A float representing the scale factor to be used for calculating the marginal changes in NpVaR.
- alpha: A float representing the significance level for calculating the NpVaR.

### Example:

```python
import numpy as np
from rmpy.NpVaR import marg_NpVaRs_scale_factor
returns = np.random.uniform(-0.05, 0.05, size=(10, 5))  # Replace this with actual returns
positions = [1000, -5000, 2000, -1000, 1500]
scale_factor = 0.2
alpha = 0.01
NpVaR_changes = marg_NpVaRs_scale_factor(returns, positions, alpha=alpha)
print(NpVaR_changes)
```
In this example, we have calculated the marginal changes in the non-parametric Aalue at Risk (NpVaR) for a portfolio consisting of 5 assets. The function takes in historical returns data for the assets, the current positions of the assets, a scale factor, and a significance level (alpha). The result is a list of the marginal changes in NpVaR for each asset, which can be useful for understanding the risk contributions of individual positions within the portfolio and designing risk management strategies to mitigate potential losses.

# 2. PaVaR Module


This module calculates Portfolio Value-at-Risk (pVaR) and associated functions for single assets and multi-asset portfolios using historical returns. It provides different Value-at-Risk (VaR) calculations, including marginal, component, and relative component VaRs. 
The functions with the `yf_` prefix use data from the Yahoo Finance API, whereas the others use the provided returns. The `_single` functions are for individual assets, and the `_port` functions are for portfolios. The `_conflevel_` functions calculate the pVaR with a specified confidence level. The `_und_vs_` functions calculate the undiversified pVaR, while the `_marginal_VaRs`, `_componets_VaRs`, and `_relcomponets_VaRs` functions calculate the marginal, component, and relative component VaRs, respectively.

Each function has its own set of input parameters, such as asset tickers, positions, start and end dates, frequency, and confidence levels. Please refer to the code above for the syntax and examples on how to use each function. The input parameters and their usage are described within the comments.

The following functions are included:

`1. yf_pVaR_single`

-  yf_pVaR_single(ticker, position, start_date, end_date,interval = 1, freq="daily", alpha=0.01, display=True)

The function 'yf_pVaR_single' enables the calculation of parametric VaR for a SINGLE POSITION by utilizing data obtained from Yahoo Finance (yFinance).

Parametric VaR is a risk management technique used to estimate the potential losses in an asset of assets, based on assumptions about the underlying distribution of the returns, such as the normal distribution.

#### Args:

- ticker: The stock symbol or identifier for the financial instrument in question (e.g. "AAPL" for Apple Inc.).
- position: The value held. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly'' " and "interval = 1", the function computes 1-month VaR).
- freq: The frequency at which returns will be downloaded.
alpha: The significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- display: A boolean value that determines whether the function should display the results of the VaR calculation. By default, this is set to True, which means that the results will be displayed.

### Example:

```python
from rmpy.PaVaR import yf_pVaR_single
ticker = ['AAPL']
position = [-1000]
start_date = '2020-01-01'
end_date = '2021-12-31'
interval = 1
freq = "daily"
alpha = 0.01
display = True
yf_pVaR_single(ticker, position, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=display)
# OR
VaR = yf_pVaR_single(ticker, position, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=False)
print(VaR)
```
This example calculates the daily parametric VaR for a short position of 1000 in Apple Inc, with a confidence level of 99%.

`2. yf_pVaR_port`

- yf_pVaR_port(tickers, positions,  start_date, end_date,interval = 1, freq="daily", alpha=0.01, display=True)

The function 'yf_pVaR_port' enables the calculation of parametric VaR for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

Parametric VaR is a statistical method of calculating the minimum amount of loss that a portfolio is likely to experience at a given confidence level, using a parametric model to estimate the distribution of returns.

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- freq: The frequency at which returns will be downloaded.
- alpha: The significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- display: A boolean value that determines whether the function should display the results of the VaR calculation. By default, this is set to True, which means that the results will be displayed.

### Example:

```python
from rmpy.PaVaR import yf_pVaR_port

tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
interval = 5
freq = "daily"
alpha = 0.01
display = True
yf_pVaR_port(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=display)
# OR
VaR = yf_pVaR_port(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=False)
print(VaR)
```
This example calculates the 5-day parametric VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

`3. pVaR_single`

- VaR = pVaR_single(returns, position, interval = 1, alpha=0.01)

The function 'pVaR_single' enables the calculation of parametric VaR for a SINGLE POSITION based on a set of returns (you can use every type of assets e.g stock, options, bonds, ecc.).

Parametric VaR is a way of estimating the potential downside risk of an asset that relies on specific assumptions about the shape and parameters of the underlying probability distribution, such as mean and standard deviation.

#### Args:

- returns: A pandas Series or NumPy array containing the historical returns of the asset or portfolio. 
- position: The size of the position in units of the asset. 
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).

### Example:

```python
from rmpy.PaVaR import pVaR_single
import numpy as np
returns = np.array([-0.02, 0.03, -0.01, 0.015, -0.002, 0.001, 0.008, 0.002, -0.006, 0.009])  # Replace this with actual returns
position = [-1000]
interval = 1
alpha = 0.01
VaR = pVaR_single(returns, position, interval=interval, alpha=alpha)
print(VaR)
```
This example calculates the parametric VaR consisting of a single short position of 1000 for the given returns.

`4. pVaR_port`

- VaR = pVaR_port(returns, position,  interval = 1, alpha=0.01)

The function 'pVaR_port' enables the calculation of parametric VaR for a PORTOFLIO based on a set of returns (you can use every type of assets e.g stock, options, bonds, ecc.).

Parametric VaR is a risk management tool that estimate the distribution of returns and potential losses in a portfolio.

#### Args:

- returns: A pandas Series or NumPy array containing the historical returns of the portfolio.
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

### Example:

```python
from rmpy.PaVaR import pVaR_port
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 2500, 7000]
interval = 1
alpha = 0.01
VaR = pVaR_port(returns, positions, interval=interval, alpha=alpha)
print(VaR)
```
This example calculates the parametric VaR consisting of a portfolio with short positions of 1000 in the first asset, long positions of 2500 in the second asset, 
and long positions of 7000 in the third asset.

`5. yf_conflevel_pVaR_single`

-  yf_conflevel_pVaR_single(ticker, position, start_date, end_date, freq="daily",interval =1,  alpha=0.01 ,confidence_level = 0.95, display=True)

The function 'yf_conflevel_pVaR_single' enables the calculation of the confidence level of parametric VaR for a SINGLE POSITION by utilizing data obtained from Yahoo Finance (yFinance).

Parametric VaR is a way to estimate the potential loss that an asset may experience over a given period, based on the assumption of a specific probability distribution and its parameters.

#### Args:

- ticker: The stock symbol or identifier for the financial instrument in question (e.g. "AAPL" for Apple Inc.).
- position: The number of shares or units held. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- freq: The frequency at which returns will be downloaded.
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: a float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- confidence_level: a float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.
- display: a boolean or string value indicating whether or not to display the results. The default value is set to True.

### Example:

```python
from rmpy.PaVaR import yf_conflevel_pVaR_single
ticker = ['AAPL']
position = [2500]
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
interval = 1
alpha = 0.01
confidence_level = 0.95
display = True
yf_conflevel_pVaR_single(ticker, position, start_date, end_date, freq=freq, interval=interval, alpha=alpha, confidence_level=confidence_level, display=display)
# OR
VaR = yf_conflevel_pVaR_single(ticker, position, start_date, end_date, freq=freq, interval=interval, alpha=alpha, confidence_level=confidence_level, display=False)
print(VaR)
```
This example calculates parametric VaR with its lower and upper bound (95% confidence) for a position of 2500 in "AAPL".

`6. yf_conflevel_pVaR_port`

-  yf_conflevel_pVaR_port(tickers, positions, start_date, end_date, freq="daily",interval =1,  alpha=0.01 ,confidence_level = 0.95, display=True)

The function 'yf_conflevel_pVaR_port' enables the calculation of the confidence level of parametric VaR for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

Parametric VaR is a measure of the worst-case loss that a portfolio is likely to experience at a given confidence level, using a parametric model to estimate the distribution of returns and potential losses.

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- freq: The frequency at which returns will be downloaded.
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- confidence_level: A float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.
- display: A boolean or string value indicating whether or not to display the results. The default value is set to True.

### Example:

```python
from rmpy.PaVaR import yf_conflevel_pVaR_port
tickers = ['AAPL', 'MSFT']
positions = [2500, -1000]
start_date = '2020-01-01'
end_date = '2021-12-31'
freq = "daily"
interval = 1
alpha = 0.01
confidence_level = 0.95
display = True
yf_conflevel_pVaR_port(tickers, positions, start_date, end_date, freq=freq, interval=interval, alpha=alpha, confidence_level=confidence_level, display=display)
# OR
VaR = yf_conflevel_pVaR_port(tickers, positions, start_date, end_date, freq=freq, interval=interval, alpha=alpha, confidence_level=confidence_level, display=False)
print(VaR)
```
This example calculates parametric VaR with its lower and upper bound (95% confidence) for a Portfolio of a position of 2500 (long position) in "AAPL" and -1000 (short position) in "MSFT.

`7. conflevel_pVaR_single`

- VaR = conflevel_pVaR_single(returns, position,interval =1,  alpha=0.01 ,confidence_level = 0.95)

The function 'conflevel_pVaR_single' enables the calculation of the confidence level of parametric VaR for a SINGLE POSITION based on a set of returns (you can use every type of assets e.g stock, options, bonds, ecc.).

Parametric VaR is a way of estimating the potential downside risk of an asset that is based on specific parametric assumptions about the underlying probability distribution, such as skewness and kurtosis.

#### Args:

- returns: A pandas Series or NumPy array containing the historical returns of the asset.
- position: The size of the position in the asset. This can be a positive or negatie number, depending if you are long or short.
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- confidence_level: A float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.
- alpha: The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

### Example:

```python
from rmpy.PaVaR import conflevel_pVaR_single
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=1000) # Replace this with your actual returns
position = [-1000]
interval = 1
confidence_level = 0.95
alpha = 0.01
VaR = conflevel_pVaR_single(returns, position, interval=interval, confidence_level=confidence_level, alpha=alpha)
print(VaR)
```
This example calculates the parametric VaR and its lower and upper bound (95% confidence interval) consisting of a position of -1000 (short position) in the given asset.

`8. conflevel_pVaR_port`

- VaR = conflevel_pVaR_port(returns, positions ,interval =1,  alpha=0.01 ,confidence_level = 0.95)

The function 'conflevel_pVaR_port' enables the calculation of the confidence level of parametric VaR for a PORTFOLIO based on a set of returns (you can use every type of assets e.g stock, options, bonds, ecc.).

Parametric VaR is a statistical method of calculating the minimum amount of loss that a portfolio is likely to experience at a given confidence level, using a parametric model to estimate the distribution of returns.

#### Args:

- returns: A pandas Series or NumPy array containing the historical returns of the portfolio.
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- confidence_level: A float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.
- alpha: The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

### Example:

```python
from rmpy.PaVaR import conflevel_pVaR_port
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 5000, -1500]
interval = 1
confidence_level = 0.95
alpha = 0.01
VaR = conflevel_pVaR_port(returns, positions, interval=interval, confidence_level=confidence_level, alpha=alpha)
print(VaR)
```
This example calculates the parametric VaR for a Portfolio and its lower and upper bound (95% confidence interval) consisting of a position of -1000 in the first asset, 5000 in the second and -1500 in the third one. 

`9. yf_und_vs_pVaR_port`

-  yf_und_vs_pVaR_port(tickers, positions, start_date, end_date,interval = 1, freq="daily", alpha=0.01, display=True)

The function 'yf_und_vs_pVaR_port' enables the calculation of the undiversified VaR and the parametric VaR for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

Undiversified VaR is a measure of the potential loss of a portfolio that does not take into account the effects of diversification, i.e., it assumes that all the assets in the portfolio move in the same direction. Parametric VaR takes into account the correlation between the assets in the portfolio, which can result in a lower estimate of potential losses. 

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- freq: The frequency at which returns will be downloaded.
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- display : A boolean or string value indicating whether or not to display the results. The default value is set to True.

### Example:

```python
from rmpy.PaVaR import yf_und_vs_pVaR_port
tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
interval = 5
freq = "daily"
alpha = 0.01
display = True
yf_und_vs_pVaR_port(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=display)
# OR
VaR = yf_und_vs_pVaR_port(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=False)
print(VaR)
```
This example calculates the 5-day undiversified and parametric VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

`10. und_vs_pVaR_port`

-  VaR = und_vs_pVaR_port(returns, position, interval = 1, alpha=0.01)

The function 'und_vs_pVaR_port' enables the calculation of the undiversified VaR and the parametric VaR for a PORTFOLIO based on a set of returns (you can use every type of assets e.g stock, options, bonds, ecc.).

Undiversified VaR is a measure of the potential loss of a portfolio that does not take into account the effects of diversification, i.e., it assumes that all the assets in the portfolio move in the same direction. Parametric VaR takes into account the correlation between the assets in the portfolio, which can result in a lower estimate of potential losses. 

#### Args:

- returns: A list or array of historical returns for a portfolio.
- positions: A list or array of current positions for the assets in the portfolio. This can be a positive or negative number, depending if you are long or short.
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

### Example:

```python
from rmpy.PaVaR import und_vs_pVaR_port
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 2500, 7000]
interval = 1
alpha = 0.01
VaR = und_vs_pVaR_port(returns, positions, interval=interval, alpha=alpha)
print(VaR)
```
This example calculates the 5-day undiversified and parametric VaR for a Portfolio with short position of 1000 in the first asset, 2500 in the second one and 7000 in the third one.

`11. yf_marginal_VaRs`

-  yf_marginal_VaRs(tickers, positions, start_date, end_date,interval = 1, freq="daily", alpha=0.01, display=display)

The function 'yf_marginal_VaRs' enables the calculation of the Marginal VaRs for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

Marginal VaR is the amount of potential loss that adding an additional unit of an asset to a portfolio would contribute to the overall VaR of the portfolio.

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- freq: The frequency at which returns will be downloaded.
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- display": A boolean or string value indicating whether or not to display the results. The default value is set to True.

### Example:

```python
from rmpy.PaVaR import yf_marginal_VaR
tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
interval = 5
freq = "daily"
alpha = 0.01
display = True
yf_marginal_VaRs(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=display)
# OR
mVars = yf_marginal_VaRs(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=False)
print(mVars)
```
This example calculates the 5-day marginal parametrics VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

`12. marginal_VaRs`

-  mVars = marginal_VaRs(returns, positions, interval=1, alpha=0.01)

The function 'marginal_VaRs' enables the calculation of the marginal VaRs of a portfolio, given a set of returns and positions.

Marginal VaR is the change in the overall VaR of a portfolio when the position size of a particular asset or position is increased or decreased.

#### Args:

- returns: A pandas Series or NumPy array containing the historical returns of the portfolio.
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

### Example:

```python
from rmpy.PaVaR import marginal_VaRs
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 2500, 7000]
interval = 1
alpha = 0.01
mVars = marginal_VaRs(returns, positions, interval=interval, alpha=alpha)
print(mVars)
```
This example calculates the marginal VaRs of a Portfolio consisting of a short positions of 1000 in the first asset, long positions of 2500 in the second asset, 
and long positions of 7000 in the third asset.

`13. yf_components_VaRs`

-  yf_components_VaRs(tickers, positions, start_date, end_date,interval = 1, freq="daily", alpha=0.01, display=True)

The function 'yf_components_VaRs' enables the calculation of the Component VaRs for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

Component VaR is a technique that breaks down the total VaR of a portfolio into the individual VaR contributions of each asset or factor, allowing investors to identify which components of the portfolio are driving the risk.

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- freq: The frequency at which returns will be downloaded.
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- display: A boolean or string value indicating whether or not to display the results. The default value is set to True.

### Example:

```python
from rmpy.PaVaR import yf_components_VaRs
tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
interval = 5
freq = "daily"
alpha = 0.01
display = True
yf_components_VaRs(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=display)
# OR
cVars = yf_components_VaRs(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=False)
print(cVars)
```
This example calculates the 5-day component parametrics VaRs for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

`14. components_VaRs`

-  cVars = components_VaRs(returns, positions, interval=1, alpha=0.01)

The function 'components_VaRs' enables the calculation of the components VaRs of a portfolio, given a set of returns.

Component VaR is a risk management tool that allows investors to analyze the contribution of each asset or factor to the overall VaR of a portfolio.

#### Args:

- returns: A pandas Series or NumPy array containing the historical returns of the portfolio.
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

### Example:

```python
from rmpy.PaVaR import components_VaRs
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 2500, 7000]
interval = 1
alpha = 0.01
cVars = components_VaRs(returns, positions, interval=interval, alpha=alpha)
print(cVars)
```
This example calculates the component VaRs of a Portfolio consisting of a short positions of 1000 in the first asset, long positions of 2500 in the second asset,  and long positions of 7000 in the third asset.

`15. yf_relcomponents_VaRs`

- yf_relcomponents_VaRs(tickers, positions, start_date, end_date,interval = 1, freq="daily", alpha=0.01, display=True)

The function 'yf_relative_components_VaRs' enables the calculation of the Relative Component VaRs for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

Relative Component VaR is a variation of Component VaR that takes into account the weight or investment in each asset or factor, allowing investors to assess the contribution of each component to the overall VaR in proportion to its size.

#### Args:

- tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- start_date: A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".
- end_date: A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".
- freq: The frequency at which returns will be downloaded.
- interval: The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.
- display: A boolean or string value indicating whether or not to display the results. The default value is set to True.

### Example:

```python
from rmpy.PaVaR import yf_relcomponents_VaRs
tickers = ['AAPL', 'MSFT']
positions = [-1000, 5000]
start_date = '2020-01-01'
end_date = '2021-12-31'
interval = 5
freq = "daily"
alpha = 0.01
display = True
yf_relcomponents_VaRs(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=display)
# OR
rcVars = yf_relcomponents_VaRs(tickers, positions, start_date, end_date, interval=interval, freq=freq, alpha=alpha, display=False)
print(rcVars)
```
This example calculates the 5-day relative component parametrics VaRs for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

`16. relcomponents_VaRs`

-  rcVars = relcomponents_VaRs(returns, positions, interval=1, alpha=0.01)

The function 'relcomponents_VaRs' enables the calculation of the relative components VaRs of a portfolio, given its returns.

Relative Component VaR is a risk management approach that breaks down the total VaR of a portfolio into the individual VaR contributions of each asset or factor, taking into account the investment or weight of each component, and provides a relative measure of the contribution of each component to the overall VaR.

#### Args:

- returns: A pandas Series or NumPy array containing the historical returns of the portfolio.
- positions: A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.
- interval: The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).
- alpha: The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.


### Example:

```python
from rmpy.PaVaR import relcomponents_VaRs
import numpy as np
returns = np.random.uniform(-0.05, 0.05, size=(10, 3))  # Replace this with actual returns
positions = [-1000, 2500, 7000]
interval = 1
alpha = 0.01
rcVars = relcomponents_VaRs(returns, positions, interval=interval, alpha=alpha)
print(rcVars)
```
This example calculates the relative component VaRs of a Portfolio consisting of a short positions of 1000 in the first asset, long positions of 2500 in the second asset, and long positions of 7000 in the third asset.

