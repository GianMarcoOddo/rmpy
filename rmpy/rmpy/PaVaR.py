
from typing import Union, List
import numpy as np

'''#################################################################### SINGLE AND PORT PARAMETRIC VAR #####################################################################'''

## 1. yf_pVaR_single #########################################################################################################################################################
############################################################################################################################################################################
 
def yf_pVaR_single(ticker : str,
                    position : Union[int,float],
                    start_date : str, 
                    end_date: str = "today",
                    interval: int = 1,
                    freq: str = "daily",
                    alpha: float = 0.01,
                    display: bool = True) -> Union[float, None]:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 
    
    Args
    -------
    The function 'yf_pVaR_single' enables the calculation of parametric VaR for a SINGLE POSITION by utilizing data obtained from Yahoo Finance (yFinance).
    
    The arguments of 'yf_pVaR_single' are:

    "ticker": The stock symbol or identifier for the financial instrument in question (e.g. "AAPL" for Apple Inc.).

    "position": The value held. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "freq": The frequency at which returns will be downloaded.

    "alpha": The significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "display": A boolean value that determines whether the function should display the results of the VaR calculation. By default, this is set to True, which means that the results will be displayed.
    
    *******************************************************************************************************************************************************************************************
    
    Example
    -------
    >>> from rmpy.PaVaR import yf_pVaR_single
    >>> ticker = ['AAPL'] 
    >>> position = [-1000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> interval = 1
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_pVaR_single(ticker, position, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=display)
    >>> # OR
    >>> VaR = yf_pVaR_single(ticker, position, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=False)
    >>> print(VaR)

    This example calculates the daily parametric VaR for a short position of 1000 in Apple Inc, with a confidence level of 99%.

    *******************************************************************************************************************************************************************************************
    '''

    import numpy as np ; import pandas as pd ; import yfinance as yf 
    from rmpy.inputs_handler import check_ticker_single, check_position_single , check_and_convert_dates , check_interval, check_freq , check_alpha , check_display

    # check ticker
    ticker = check_ticker_single(ticker, (str, list, np.ndarray, pd.Series), "ticker", single_value=True)
    # check position
    position = check_position_single(position)
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS

    def get_data(ticker, start_date, end_date, freq):

        freq_to_interval = {"daily": "1d","weekly": "1wk","monthly": "1mo"}
        interval = freq_to_interval.get(freq)
        try:
            data = yf.download(ticker, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        except:
            raise Exception("The download process encountered an error. Please verify your internet connection")
        if data.empty or data.shape[0] < 3:
            raise Exception(f"Either the provided ticker {ticker} may have been delisted or the date range is too small for the VaR calculation")
        
        # Calculate returns and drop missing values
        returns = data["Adj Close"].pct_change().dropna()

        return returns

    # Call function with your parameters
    returns = get_data(ticker, start_date, end_date, freq)

    VaR = pVaR_single(returns, position, interval = interval,  alpha = alpha)

    # Print settings
    conf_level = int((1-alpha)*100)
    display = str(display).lower() in ['true', 't']
    if display:
            print(f"Parametric VaR({conf_level}%,{interval}) for {ticker} =  {max(0,VaR):.3f} $")
    return VaR

## 2. yf_pVaR_port #########################################################################################################################################################
########################################################################################################################################################################

def yf_pVaR_port(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  interval: int = 1,
                  freq: str = "daily",
                  alpha: float = 0.01,
                  display: bool = True) -> Union[float, None]:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    -------
    The function 'yf_pVaR_port' enables the calculation of parametric VaR for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

    The arguments of 'yf_pVaR_port' are:

    "tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "freq": The frequency at which returns will be downloaded.

    "alpha": The significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "display": A boolean value that determines whether the function should display the results of the VaR calculation. By default, this is set to True, which means that the results will be displayed.

    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import yf_pVaR_port
    >>> tickers = ['AAPL','MSFT'] 
    >>> positions = [-1000,5000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> interval = 5
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_pVaR_port(tickers, positions, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=display)
    >>> # OR
    >>> VaR = yf_pVaR_port(ticker, position, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=False)
    >>> print(VaR)

    This example calculates the 5-day parametric VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

    *******************************************************************************************************************************************************************************************
    '''

    import yfinance as yf
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_and_convert_dates , check_interval, check_freq , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS  
    
    def get_data(tickers, start_date, end_date, freq):
        # Create a dictionary for frequency-to-interval mapping
        freq_to_interval = {"daily": "1d", "weekly": "1wk", "monthly": "1mo"}
        # Perform the download once using the relevant interval
        interval = freq_to_interval[freq]
        data = yf.download(tickers, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        if data.empty or data.shape[0] < 3:
            raise Exception(f"The date range is too small for the VaR calculation or all the provided tickers are wrong")
        # Process the downloaded data
        data = data["Adj Close"].dropna(axis=1, how='all')
        returns = data.pct_change().dropna()
        # Handle failed tickers
        ticker_downloaded = returns.columns
        failed_tickers = [ticker for ticker in tickers if ticker not in ticker_downloaded]
        if failed_tickers:
            tickers_str = ', '.join(failed_tickers)
            is_are = 'is' if len(failed_tickers) == 1 else 'are'
            raise Exception(f"{tickers_str} {is_are} incorrect or may have been delisted. Please check if these tickers are correct")
        return returns
    
    returns = get_data(tickers, start_date, end_date, freq)

    from rmpy.PaVaR import pVaR_port

    VaR = pVaR_port(returns, positions, interval = interval,  alpha = alpha)

    # Print settings
    conf_level = int((1-alpha)*100)
    display = str(display).lower() in ['true', 't']
    if display:
            print(f"Parametric VaR({conf_level}%,{interval}) for the Portfolio =  {max(0,VaR):.3f} $")
    return VaR


## 3. pVaR_single #########################################################################################################################################################
########################################################################################################################################################################


def pVaR_single(returns : np.ndarray,
                position: Union[int,float],
                interval: int = 1,
                alpha: float = 0.01) -> Union[float, None]:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    -------
    The function 'pVaR_single' enables the calculation of parametric VaR for a SINGLE POSITION  based on a set of returns.

    The arguments of the function are: 

    "returns": a pandas Series or NumPy array containing the historical returns of the asset or portfolio. 

    "position": the size of the position in units of the asset. 

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. 
    By default, this is set to 0.01 (1%).

    VaR is then calculated as the product of the z-score q (depending on the alpha parameter), the position (investment size), the standard deviation of returns, 
    and the square root of the interval (time horizon). 
    The calculation assumes that the returns are normally distributed, and the position's risk can be captured by the standard deviation of the returns.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import pVaR_single ; import numpy as np
    >>> returns = np.array([-0.02, 0.03, -0.01, 0.015, -0.002, 0.001, 0.008, 0.002, -0.006, 0.009]) # Replace this with actual returns
    >>> position = [-1000]
    >>> interval = 1
    >>> alpha = 0.01
    >>> VaR = pVaR_single(returns, position, interval = interval, alpha=alpha)
    >>> print(VaR)

    This example calculates the parametric VaR consisting of a single short position of 1000 for the given returns.

    *****************************************************************************************************************************************************************************************************************
    '''

    import pandas as pd ; import numpy as np ; from scipy.stats import norm
    from rmpy.inputs_handler import validate_returns_single , check_position_single, check_interval, check_alpha

    # check returns
    returns = validate_returns_single(returns)
    # check position
    position = check_position_single(position) 
    if position < 0: 
        position = -1 * position
    # check interval 
    interval = check_interval(interval)
    # check alpha
    alpha = check_alpha(alpha)

    sigma = np.std(returns)
    q = norm.ppf(1-alpha, loc=0, scale=1)
    VaR = q * position * sigma * np.sqrt(interval)

    return max(0,VaR)

## 4. pVaR_port #########################################################################################################################################################
########################################################################################################################################################################

def pVaR_port(returns: np.ndarray,
            positions: List[Union[int, float]],
            interval: int = 1,
            alpha: float = 0.01)-> Union[float, None]:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 
    
    Args
    -------
    The function 'pVaR_port' enables the calculation of parametric VaR for a PORTOFLIO based on a set of returns.

    The arguments of the function are: 

    "returns": a pandas Series or NumPy array containing the historical returns of the portfolio.

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

    VaR: Calculates the Value at Risk for the portfolio as the product of the z-score q(depending on the alpha parameter) , the portfolio standard deviation sigma (Variance - Covariance Matrix), 
    and the square root of the interval (time horizon). 
    The calculation assumes that the portfolio's returns are normally distributed and the portfolio risk can be captured by the standard deviation of the returns.
    The function returns the maximum of 0 or the calculated VaR. This ensures that the Value at Risk is a non-negative number, as VaR is meant to represent the potential loss.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import pVaR_port ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [-1000,2500,7000]
    >>> interval = 1
    >>> alpha = 0.01
    >>> VaR = pVaR_port(returns, positions, interval = interval, alpha=alpha)

    This example calculates the parametric VaR consisting of a portfolio with short positions of 1000 in the first asset, long positions of 2500 in the second asset, 
    and long positions of 7000 in the third asset.
    
    ************************************************************************************************************************************************************************************************************************
    '''

    import numpy as np ; from scipy.stats import norm
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_interval, check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check interval 
    interval = check_interval(interval)
    # check alpha
    alpha = check_alpha(alpha)
    
    variance_cov_matrix = np.cov(returns, rowvar=False)
    var_port = np.dot(np.dot(positions,variance_cov_matrix),positions.T) 
    sigma = np.sqrt(var_port)
    q = norm.ppf(1-alpha, loc=0, scale=1)
    VaR = q * sigma * np.sqrt(interval)

    return max(0,VaR)

'''########################################################## CONFIDENCE LEVEL SINGLE AND PORT PARAMETRIC VAR #############################################################'''

## 5. yf_conflevel_pVaR_single ###############################################################################################################################################
############################################################################################################################################################################

def yf_conflevel_pVaR_single(ticker : str,
                    position : Union[int,float],
                    start_date : str, 
                    end_date: str = "today",
                    freq: str = "daily",
                    interval : int = 1,
                    alpha: float = 0.01,
                    confidence_level: float = 0.95,
                    display: bool = True) -> dict:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    -------
    The function 'yf_conflevel_pVaR_single' enables the calculation of the confidence level of parametric VaR for a SINGLE POSITION by utilizing data obtained from Yahoo Finance (yFinance).

    The arguments of the function are:

    "ticker": The stock symbol or identifier for the financial instrument in question (e.g. "AAPL" for Apple Inc.).

    "position": The number of shares or units held. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "freq": The frequency at which returns will be downloaded.

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "confidence_level": A float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.

    "display": A boolean or string value indicating whether or not to display the results. The default value is set to True.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import yf_conflevel_pVaR_single
    >>> ticker = ['AAPL']
    >>> position = [2500]
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> interval = 1
    >>> alpha = 0.01
    >>> confidence_level = 0.95
    >>> display = True
    >>> yf_conflevel_pVaR_single(ticker, position, start_date, end_date, freq=freq,interval =interval,  alpha=alpha ,confidence_level = confidence_level, display=display)
    >>> # OR
    >>> VaR = yf_conflevel_pVaR_single(ticker, position, start_date, end_date, freq=freq,interval =interval, alpha=alpha ,confidence_level = confidence_level, display=False)
    >>> print(VaR)

    This example calculates parametric VaR with its lower and upper bound (95% confidence) for a position of 2500 in "AAPL".

    *********************************************************************************************************************************************************************************************************
    '''

    import numpy as np ; import pandas as pd ; import yfinance as yf
    from rmpy.inputs_handler import check_ticker_single, check_position_single , check_and_convert_dates , check_freq ,check_interval , check_alpha ,check_confidence_level, check_display

    # check ticker
    ticker = check_ticker_single(ticker, (str, list, np.ndarray, pd.Series), "ticker", single_value=True)
    # check position
    position = check_position_single(position) 
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check confidence_level 
    confidence_level = check_confidence_level(confidence_level)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])
    
    # DOWNLOADING PROCESS

    def get_data(ticker, start_date, end_date, freq):

        freq_to_interval = {"daily": "1d","weekly": "1wk","monthly": "1mo"}
        interval = freq_to_interval.get(freq)
        try:
            data = yf.download(ticker, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        except:
            raise Exception("The download process encountered an error. Please verify your internet connection")
        if data.empty:
            raise Exception(f"The provided ticker {ticker} may have been delisted or it might be wrong")

        # Calculate returns and drop missing values
        returns = data["Adj Close"].pct_change().dropna()

        # Raise an exception if less than two returns have been downloaded
        if returns.shape[0] < 2:
            raise Exception(f"Insufficient data to calculate the Parametric VaR for {ticker}")
        return returns

    # Call function with your parameters
    returns = get_data(ticker, start_date, end_date, freq)

    from rmpy.PaVaR import conflevel_pVaR_single
    summary = conflevel_pVaR_single(returns, position, interval = interval, confidence_level = confidence_level, alpha = alpha)
    # Print settings
    display = str(display).lower() in ['true', 't']
    if display:
            rounded_summary = {key: round(value, 4) for key, value in summary.items()}
            print(rounded_summary)
    return summary


## 6. yf_conflevel_pVaR_port ###############################################################################################################################################
############################################################################################################################################################################

def yf_conflevel_pVaR_port(tickers : List[str],
                    positions : List[Union[int,float]],
                    start_date : str, 
                    end_date: str = "today",
                    freq: str = "daily",
                    interval : int = 1,
                    alpha: float = 0.01,
                    confidence_level: float = 0.95,
                    display: bool = True) -> dict:
    
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    -------
    The function 'yf_conflevel_pVaR_port' enables the calculation of the confidence level of parametric VaR for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

    The arguments of the function are: 

    "tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "freq": The frequency at which returns will be downloaded.

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "confidence_level": A float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.

    "display": A boolean or string value indicating whether or not to display the results. The default value is set to True.

    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import yf_conflevel_pVaR_port
    >>> tickers = ['AAPL',"MSFT"]
    >>> positions = [2500,-1000]
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> interval = 1
    >>> alpha = 0.01
    >>> confidence_level = 0.95
    >>> display = True
    >>> yf_conflevel_pVaR_port(tickers, positions, start_date, end_date, freq=freq,interval =interval,  alpha=alpha ,confidence_level = confidence_level, display=display)
    >>> # OR
    >>> VaR = yf_conflevel_pVaR_port(tickers, positions, start_date, end_date, freq=freq,interval =interval, alpha=alpha ,confidence_level = confidence_level, display=False)
    >>> print(VaR)

    This example calculates parametric VaR with its lower and upper bound (95% confidence) for a Portfolio of a position of 2500 in "AAPL" and - 1000 in "MSFT.

    *******************************************************************************************************************************************************************************************
    '''
   
    import yfinance as yf
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_and_convert_dates , check_interval, check_freq , check_alpha ,check_confidence_level, check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check confidence_level 
    confidence_level = check_confidence_level(confidence_level)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS  
    
    def get_data(tickers, start_date, end_date, freq):
        # Create a dictionary for frequency-to-interval mapping
        freq_to_interval = {"daily": "1d", "weekly": "1wk", "monthly": "1mo"}
        # Perform the download once using the relevant interval
        interval = freq_to_interval[freq]
        data = yf.download(tickers, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        if data.empty or data.shape[0] < 3:
            raise Exception(f"The date range is too small for the VaR calculation")
        # Process the downloaded data
        data = data["Adj Close"].dropna(axis=1, how='all')
        returns = data.pct_change().dropna()
        # Handle failed tickers
        ticker_downloaded = returns.columns
        failed_tickers = [ticker for ticker in tickers if ticker not in ticker_downloaded]
        if failed_tickers:
            tickers_str = ', '.join(failed_tickers)
            is_are = 'is' if len(failed_tickers) == 1 else 'are'
            raise Exception(f"{tickers_str} {is_are} incorrect or may have been delisted. Please check if these tickers are correct")
        return returns
    
    returns = get_data(tickers, start_date, end_date, freq)

    from rmpy.PaVaR import conflevel_pVaR_port
    summary = conflevel_pVaR_port(returns, positions, interval = interval,confidence_level= confidence_level,  alpha = alpha)
    # Print settings
    display = str(display).lower() in ['true', 't']
    if display:
            rounded_summary = {key: round(value, 4) for key, value in summary.items()}
            print(rounded_summary)
    return summary

## 7. conflevel_pVaR_single ###############################################################################################################################################
############################################################################################################################################################################

def conflevel_pVaR_single(returns: np.ndarray,
                        position = List[Union[int,float]],
                        interval: int = 1,
                        confidence_level: float = 0.95,
                        alpha: float = 0.01) -> dict:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'conflevel_pVaR_single' enables the calculation of the confidence level of parametric VaR for a SINGLE POSITION based on a set of returns.

    The arguments of the function are: 

    "returns": a pandas Series or NumPy array containing the historical returns of the asset.

    "position": The size of the position in the asset. This can be a positive or negative number, depending if you are long or short.

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "confidence_level": A float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.

    "alpha": The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import conflevel_pVaR_single ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=1000)  # Replace this with your actual returns
    >>> position = [-1000]
    >>> interval = 1
    >>> confidence_level = 0.95
    >>> alpha = 0.01
    >>> VaR = conflevel_pVaR_single(returns, position,interval = interval, confidence_level= confidence_level, alpha=alpha)
    >>> print(VaR)

    This example calculates the parametric VaR and its lower and upper bound (95% confidence interval) consisting of a position of -1000 in the given asset.

    *******************************************************************************************************************************************************************************************
    '''

    import pandas as pd ; import numpy as np ; from scipy.stats import norm
    from rmpy.inputs_handler import validate_returns_single , check_position_single, check_interval, check_confidence_level, check_alpha

    # check returns
    returns = validate_returns_single(returns)
    # check position
    position = check_position_single(position) 
    if position <0:
        position = -1* position
    # check interval
    interval = check_interval(interval)
    # check confidence_level 
    confidence_level = check_confidence_level(confidence_level)
    # check alpha
    alpha = check_alpha(alpha)

    from rmpy.PaVaR import pVaR_single

    VaR = pVaR_single(returns, position, interval =interval,  alpha = alpha)

    sigma = np.std(returns)

    q = norm.ppf(1-alpha, loc=0, scale=1)
    se = sigma / (np.sqrt(returns.shape[0]))
    lower_bound = (q * position * (sigma - 2 * se)) * np.sqrt(interval)
    upper_bound = (q * position * (sigma + 2 * se)) * np.sqrt(interval)

    summary = {"pVaR":VaR,"lower_bound":lower_bound,"upper_bound":upper_bound}
    return summary


## 8. conflevel_pVaR_port ###############################################################################################################################################
############################################################################################################################################################################

def conflevel_pVaR_port(returns: np.ndarray,
                        positions: List[Union[int,str]],
                        interval: int = 1,
                        confidence_level: float = 0.95 ,
                        alpha: float = 0.01) -> dict:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'conflevel_pVaR_port' enables the calculation of the confidence level of parametric VaR for a PORTFOLIO based on a set of returns.

    The arguments of the function are: 

    "returns": A pandas Series or NumPy array containing the historical returns of the portfolio.

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "confidence_level": A float specifying the confidence level for the VaR calculation. By default, it is set to 0.95.

    "alpha": The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import conflevel_pVaR_port ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> position = [-1000,5000,-1500]
    >>> interval = 1
    >>> confidence_level = 0.95
    >>> alpha = 0.01
    >>> VaR = conflevel_pVaR_port(returns, positions ,interval = interval, confidence_level= confidence_level, alpha=alpha)
    >>> print(VaR)

    This example calculates the parametric VaR for a Portfolio and its lower and upper bound (95% confidence interval) consisting of a position of -1000 in the first asset, 5000 in the second and
    -1500 in the third one. 

    *******************************************************************************************************************************************************************************************
    '''

    import numpy as np ; from scipy.stats import norm
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_interval, check_confidence_level,  check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check interval
    interval = check_interval(interval)
    # check confidence_level 
    confidence_level = check_confidence_level(confidence_level)
    # check alpha
    alpha = check_alpha(alpha)
    
    from rmpy.PaVaR import pVaR_port

    VaR = pVaR_port(returns, positions, interval = interval,  alpha = alpha)

    variance_cov_matrix = np.cov(returns, rowvar=False)
    var_port = np.dot(np.dot(positions,variance_cov_matrix),positions.T) 
    sigma = np.sqrt(var_port)

    q = norm.ppf(1-alpha, loc=0, scale=1)
    se = sigma / (np.sqrt(returns.shape[0]))
    lower_bound = (q * (sigma - 2 * se)) * np.sqrt(interval)
    upper_bound = (q * (sigma + 2 * se)) * np.sqrt(interval)

    summary = {"pVaR":VaR,"lower_bound":lower_bound,"upper_bound":upper_bound}
    return summary


'''########################################################## UNDIVERSIFIED VAR VS PORT VAR #############################################################'''


## 9. yf_und_vs_pVaR_port ###############################################################################################################################################
##################################################################################################################################################################

def yf_und_vs_pVaR_port(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  interval: int = 1,
                  freq: str = "daily",
                  alpha: float = 0.01,
                  display: bool = True) -> dict:
    """
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'yf_und_vs_pVaR_port' enables the calculation of the undiversified VaR and the parametric VaR for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

    The arguments of the function are: 

    "tickers": A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "freq": The frequency at which returns will be downloaded.

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "display": A boolean or string value indicating whether or not to display the results. The default value is set to True.

    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import yf_und_vs_pVaR_port
    >>> tickers = ['AAPL','MSFT'] 
    >>> positions = [-1000,5000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> interval = 5
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_und_vs_pVaR_port(tickers, positions, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=display)
    >>> # OR
    >>> VaR = yf_und_vs_pVaR_port(ticker, position, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=False)
    >>> print(VaR)

    This example calculates the 5-day undiversified and parametric VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

    *******************************************************************************************************************************************************************************************
    """
    
    import yfinance as yf
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_and_convert_dates , check_interval, check_freq , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS  
    
    def get_data(tickers, start_date, end_date, freq):
        # Create a dictionary for frequency-to-interval mapping
        freq_to_interval = {"daily": "1d", "weekly": "1wk", "monthly": "1mo"}
        # Perform the download once using the relevant interval
        interval = freq_to_interval[freq]
        data = yf.download(tickers, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        if data.empty or data.shape[0] < 3:
            raise Exception(f"The date range is too small for the VaR calculation")
        # Process the downloaded data
        data = data["Adj Close"].dropna(axis=1, how='all')
        returns = data.pct_change().dropna()
        # Handle failed tickers
        ticker_downloaded = returns.columns
        failed_tickers = [ticker for ticker in tickers if ticker not in ticker_downloaded]
        if failed_tickers:
            tickers_str = ', '.join(failed_tickers)
            is_are = 'is' if len(failed_tickers) == 1 else 'are'
            raise Exception(f"{tickers_str} {is_are} incorrect or may have been delisted. Please check if these tickers are correct")
        return returns
    
    returns = get_data(tickers, start_date, end_date, freq)

    from rmpy.PaVaR import und_vs_pVaR_port

    summary = und_vs_pVaR_port(returns, positions, interval = interval,  alpha = alpha)

    # Print settings
    conf_level = int((1-alpha)*100)
    display = str(display).lower() in ['true', 't']
    if display:
            rounded_summary = {key: round(value, 4) for key, value in summary.items()}
            print(rounded_summary)
    return summary

## 10. und_vs_pVaR_port ###############################################################################################################################################
##################################################################################################################################################################

def und_vs_pVaR_port(returns: np.ndarray,
                    positions: List[Union[int,float]],
                    interval: int = 1,
                    alpha: float = 0.01) -> dict:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'und_vs_pVaR_port' enables the calculation of the undiversified VaR and the parametric VaR for a PORTFOLIO based on a set of returns.

    The arguments of the function are:

    "returns": A list or array of historical returns for a portfolio.

    "positions": A list or array of current positions for the assets in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import und_vs_pVaR_port ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [-1000,2500,7000]
    >>> interval = 1
    >>> alpha = 0.01
    >>> VaR = und_vs_pVaR_port(returns, positions, interval = interval, alpha=alpha)
    >>> print(VaR)

    This example calculates the 5-day undiversified and parametric VaR for a Portfolio with short position of 1000 in the first asset, 2500 in the second one and 7000 in the third one.
    
    ************************************************************************************************************************************************************************************************************************

    '''

    import numpy as np 
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_interval, check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check interval 
    interval = check_interval(interval)
    # check alpha
    alpha = check_alpha(alpha)
    # check display

    # Undiversified VaR

    from rmpy.PaVaR import pVaR_single

    single_VaRs = []

    for index in range(0,returns.shape[1]):
        VaR = pVaR_single(returns[:, index],np.array([positions[index]]),interval = interval,  alpha = alpha)
        single_VaRs.append(VaR)

    Und_VaR = np.sum(single_VaRs)

    from rmpy.PaVaR import pVaR_port

    VaR = pVaR_port(returns,positions,interval = interval,  alpha = alpha)

    summary = {"undVaR":Und_VaR,"VaR":VaR}
    return summary

'''########################################################## MARG - COMP and RELCOMP VAR #############################################################'''


## 11. yf_marginal_VaRs ###############################################################################################################################################
##################################################################################################################################################################

def yf_marginal_VaRs(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  interval: int = 1,
                  freq: str = "daily",
                  alpha: float = 0.01,
                  display: bool = True) ->list:
    
    """
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'yf_marginal_VaRs' enables the calculation of the Marginal VaRs for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).
    
    The arguments of the function are: 

    "tickers": A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "freq": The frequency at which returns will be downloaded.

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "display": A boolean or string value indicating whether or not to display the results. The default value is set to True.

    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import yf_marginal_VaRs
    >>> tickers = ['AAPL','MSFT'] 
    >>> positions = [-1000,5000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> interval = 5
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_marginal_VaRs(tickers, positions, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=display)
    >>> # OR
    >>> mVars = yf_marginal_VaRs(ticker, position, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=False)
    >>> print(mVars)

    This example calculates the 5-day marginal parametrics VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

    *******************************************************************************************************************************************************************************************
    """
    import yfinance as yf
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_and_convert_dates , check_interval, check_freq , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS  
    
    def get_data(tickers, start_date, end_date, freq):
        # Create a dictionary for frequency-to-interval mapping
        freq_to_interval = {"daily": "1d", "weekly": "1wk", "monthly": "1mo"}
        # Perform the download once using the relevant interval
        interval = freq_to_interval[freq]
        data = yf.download(tickers, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        if data.empty or data.shape[0] < 3:
            raise Exception(f"The date range is too small for the VaR calculation")
        # Process the downloaded data
        data = data["Adj Close"].dropna(axis=1, how='all')
        returns = data.pct_change().dropna()
        # Handle failed tickers
        ticker_downloaded = returns.columns
        failed_tickers = [ticker for ticker in tickers if ticker not in ticker_downloaded]
        if failed_tickers:
            tickers_str = ', '.join(failed_tickers)
            is_are = 'is' if len(failed_tickers) == 1 else 'are'
            raise Exception(f"{tickers_str} {is_are} incorrect or may have been delisted. Please check if these tickers are correct")
        return returns
    
    returns = get_data(tickers, start_date, end_date, freq)

    from rmpy.PaVaR import marginal_VaRs

    marginal_VaRs = marginal_VaRs(returns, positions, interval = interval,  alpha = alpha)
    if display:
            rounded_marginal_VaRs = [round(value, 4) for value in marginal_VaRs]
            var_dict = dict(zip(tickers, rounded_marginal_VaRs))
            print(var_dict)
    return marginal_VaRs

## 12. marginal_VaRs ###############################################################################################################################################
##################################################################################################################################################################

def marginal_VaRs(returns: np.ndarray,
            positions: List[Union[int, float]],
            interval: int = 1,
            alpha: float = 0.01)  -> list:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'marginal_VaRs' enables the calculation of the marginal VaRs of a portfolio, given a set of returns and positions

    The arguments of the function are: 

    "returns": A pandas Series or NumPy array containing the historical returns of the portfolio.

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

    *******************************************************************************************************************************************************************************************

    Example
    ------- 
    >>> from rmpy.PaVaR import marginal_VaRs ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [-1000,2500,7000]
    >>> interval = 1
    >>> alpha = 0.01
    >>> mVars = marginal_VaRs(returns, positions, interval=interval, alpha=alpha)
    >>> print(mVars)

    This example calculates the marginal VaRs of a Portfolio consisting of a short positions of 1000 in the first asset, long positions of 2500 in the second asset, 
    and long positions of 7000 in the third asset.
    
    ************************************************************************************************************************************************************************************************************************ 
    '''

    import numpy as np ; from scipy.stats import norm
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_interval, check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check interval 
    interval = check_interval(interval)
    # check alpha
    alpha = check_alpha(alpha)
    
    variance_cov_matrix = np.cov(returns, rowvar=False)
    var_port = np.dot(np.dot(positions,variance_cov_matrix),positions.T) 
    sigma = np.sqrt(var_port)
    q = norm.ppf(1-alpha, loc=0, scale=1)
    VaR = q * sigma * np.sqrt(interval)
    marginal_VaRs =  ((np.dot(variance_cov_matrix,positions.T)) / var_port) * VaR

    return marginal_VaRs


## 13. yf_components_VaRs ###############################################################################################################################################
##################################################################################################################################################################

def yf_components_VaRs(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  interval: int = 1,
                  freq: str = "daily",
                  alpha: float = 0.01,
                  display: bool = True) ->list:
    """
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'yf_components_VaRs' enables the calculation of the Component VaRs for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

    The arguments of the function are: 

    "tickers": A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "freq": The frequency at which returns will be downloaded.

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "display": A boolean or string value indicating whether or not to display the results. The default value is set to True.
    
    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import yf_components_VaRs
    >>> tickers = ['AAPL','MSFT'] 
    >>> positions = [-1000,5000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> interval = 5
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_components_VaRs(tickers, positions, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=display)
    >>> # OR
    >>> cVars = yf_components_VaRs(ticker, position, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=False)
    >>> print(cVars)

    This example calculates the 5-day component parametrics VaRs for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

    *******************************************************************************************************************************************************************************************
    """
    import yfinance as yf
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_and_convert_dates , check_interval, check_freq , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS  
    
    def get_data(tickers, start_date, end_date, freq):
        # Create a dictionary for frequency-to-interval mapping
        freq_to_interval = {"daily": "1d", "weekly": "1wk", "monthly": "1mo"}
        # Perform the download once using the relevant interval
        interval = freq_to_interval[freq]
        data = yf.download(tickers, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        if data.empty or data.shape[0] < 3:
            raise Exception(f"The date range is too small for the VaR calculation")
        # Process the downloaded data
        data = data["Adj Close"].dropna(axis=1, how='all')
        returns = data.pct_change().dropna()
        # Handle failed tickers
        ticker_downloaded = returns.columns
        failed_tickers = [ticker for ticker in tickers if ticker not in ticker_downloaded]
        if failed_tickers:
            tickers_str = ', '.join(failed_tickers)
            is_are = 'is' if len(failed_tickers) == 1 else 'are'
            raise Exception(f"{tickers_str} {is_are} incorrect or may have been delisted. Please check if these tickers are correct")
        return returns
    
    returns = get_data(tickers, start_date, end_date, freq)

    from rmpy.PaVaR import components_VaRs

    components_VaRs = components_VaRs(returns, positions, interval = interval,  alpha = alpha)
    if display:
            rounded_components_VaRs = [round(value, 4) for value in components_VaRs]
            var_dict = dict(zip(tickers, rounded_components_VaRs))
            print(var_dict)
    return components_VaRs


## 14. components_VaRs ###############################################################################################################################################
##################################################################################################################################################################

def components_VaRs(returns: np.ndarray,
            positions: List[Union[int, float]],
            interval: int = 1,
            alpha: float = 0.01) ->list:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'components_VaRs' enables the calculation of the components VaRs of a portfolio, given a set of returns

    The arguments of the function are: 

    "returns": A pandas Series or NumPy array containing the historical returns of the portfolio.

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.

    
    *******************************************************************************************************************************************************************************************

    Example
    ------- 
    >>> from rmpy.PaVaR import components_VaRs ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [-1000,2500,7000]
    >>> interval = 1
    >>> alpha = 0.01
    >>> cVars = components_VaRs(returns, positions, interval=interval, alpha=alpha)
    >>> print(cVars)

    This example calculates the componet VaRs of a Portfolio consisting of a short positions of 1000 in the first asset, long positions of 2500 in the second asset, 
    and long positions of 7000 in the third asset.
    
    ************************************************************************************************************************************************************************************************************************ 
    '''

    import numpy as np ; from scipy.stats import norm
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_interval, check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check interval 
    interval = check_interval(interval)
    # check alpha
    alpha = check_alpha(alpha)
    
    variance_cov_matrix = np.cov(returns, rowvar=False)
    var_port = np.dot(np.dot(positions,variance_cov_matrix),positions.T) 
    sigma = np.sqrt(var_port)
    q = norm.ppf(1-alpha, loc=0, scale=1)
    VaR = q * sigma * np.sqrt(interval)
    marginal_VaRs =  ((np.dot(variance_cov_matrix,positions.T)) / var_port) * VaR
    components_VaRs = marginal_VaRs.T * positions

    return components_VaRs

## 15. yf_relcomponents_VaRs ###############################################################################################################################################
##################################################################################################################################################################

def yf_relcomponents_VaRs(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  interval: int = 1,
                  freq: str = "daily",
                  alpha: float = 0.01,
                  display: bool = True)  ->list:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'yf_relative_components_VaRs' enables the calculation of the Relative Component VaRs for a PORTFOLIO by utilizing data obtained from Yahoo Finance (yFinance).

    The arguments of the function are: 

    "tickers: A list of strings representing the tickers of the assets in the portoflio. #### Note that the tickers provided should be part of the portoflio ####

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "start_date": A string representing the starting date for the historical data used in the VaR calculation. This should be in the format "YYYY-MM-DD".

    "end_date": A string representing the ending date for the historical data used in the VaR calculation. This should also be in the format "YYYY-MM-DD".By default, this is set to "today".

    "freq": The frequency at which returns will be downloaded.

    "interval": The time horizon of VaR. It is related to frequency (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": A float specifying the significance level for the VaR calculation. By default, this is set to 0.01, which corresponds to a 99% confidence level.

    "display": A boolean or string value indicating whether or not to display the results. The default value is set to True.
    
    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.PaVaR import yf_relcomponents_VaRs
    >>> tickers = ['AAPL','MSFT'] 
    >>> positions = [-1000,5000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> interval = 5
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_relcomponents_VaRs(tickers, positions, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=display)
    >>> # OR
    >>> rcVars = yf_relcomponents_VaRs(ticker, position, start_date, end_date,interval = interval, freq=freq, alpha=alpha, display=False)
    >>> print(rcVars)

    This example calculates the 5-day relative component parametrics VaRs for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

    *******************************************************************************************************************************************************************************************
    '''
    import yfinance as yf
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_and_convert_dates , check_interval, check_freq , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check interval 
    interval = check_interval(interval)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS  
    
    def get_data(tickers, start_date, end_date, freq):
        # Create a dictionary for frequency-to-interval mapping
        freq_to_interval = {"daily": "1d", "weekly": "1wk", "monthly": "1mo"}
        # Perform the download once using the relevant interval
        interval = freq_to_interval[freq]
        data = yf.download(tickers, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        if data.empty or data.shape[0] < 3:
            raise Exception(f"The date range is too small for the VaR calculation")
        # Process the downloaded data
        data = data["Adj Close"].dropna(axis=1, how='all')
        returns = data.pct_change().dropna()
        # Handle failed tickers
        ticker_downloaded = returns.columns
        failed_tickers = [ticker for ticker in tickers if ticker not in ticker_downloaded]
        if failed_tickers:
            tickers_str = ', '.join(failed_tickers)
            is_are = 'is' if len(failed_tickers) == 1 else 'are'
            raise Exception(f"{tickers_str} {is_are} incorrect or may have been delisted. Please check if these tickers are correct")
        return returns
    
    returns = get_data(tickers, start_date, end_date, freq)

    from rmpy.PaVaR import relcomponents_VaRs

    relcomponents_VaRs = relcomponents_VaRs(returns, positions, interval = interval,  alpha = alpha)
    if display:
            rounded_relcomponents_VaRs = [round(value, 4) for value in relcomponents_VaRs]
            var_dict = dict(zip(tickers, rounded_relcomponents_VaRs))
            print(var_dict)
    return relcomponents_VaRs


## 16. relcomponents_VaRs ###############################################################################################################################################
##################################################################################################################################################################


def relcomponents_VaRs(returns: np.ndarray,
            positions: List[Union[int, float]],
            interval: int = 1,
            alpha: float = 0.01)  ->list:
    '''
    *******************************************************************************************************************************************************************************************
   
    DEFINITION
    -------
    Parametric VaR (Value at Risk) is risk management tool that measures the maximum potential loss that an investment or portfolio may experience over a given time period, with a given level of 
    confidence. Parametric VaR assumes the returns of an investment or portfolio follow a normal distribution.

    The parametric VaR approach is used in situations where the returns of an investment is assumed to follow a normal distribution. It is important to note that this approach may not
    be appropriate for investments with non-normal returns.

    ******************************************************************************************************************************************************************************************* 

    Args
    ------- 
    The function 'relcomponents_VaRs' enables the calculation of the relative components VaRs of a portfolio 

    The arguments of the function are: 

    "returns": A pandas Series or NumPy array containing the historical returns of the portfolio.

    "positions": A list of integers or floats representing each position in the portfolio. This can be a positive or negative number, depending if you are long or short.

    "interval": The time horizon of VaR. It is related to frequency of data used (e.g. if "freq = 'monthly' " and "interval = 1", the function computes 1-month VaR).

    "alpha": The level of confidence for the VaR calculation. By default is set to 0.01, which represents a 99% confidence level.
    
    *******************************************************************************************************************************************************************************************

    Example
    ------- 
    >>> from rmpy.PaVaR import relcomponents_VaRs ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [-1000,2500,7000]
    >>> interval = 1
    >>> alpha = 0.01
    >>> rcVars = relcomponents_VaRs(returns, positions, interval=interval, alpha=alpha)
    >>> print(rcVars)

    This example calculates the relative componet VaRs of a Portfolio consisting of a short positions of 1000 in the first asset, long positions of 2500 in the second asset, 
    and long positions of 7000 in the third asset.
    
    ************************************************************************************************************************************************************************************************************************ 
    '''

    import numpy as np ; from scipy.stats import norm
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_interval, check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check interval 
    interval = check_interval(interval)
    # check alpha
    alpha = check_alpha(alpha)
    
    variance_cov_matrix = np.cov(returns, rowvar=False)
    var_port = np.dot(np.dot(positions,variance_cov_matrix),positions.T) 
    sigma = np.sqrt(var_port)
    q = norm.ppf(1-alpha, loc=0, scale=1)
    VaR = q * sigma * np.sqrt(interval)
    marginal_VaRs =  ((np.dot(variance_cov_matrix,positions.T)) / var_port) * VaR
    components_VaRs = marginal_VaRs.T * positions
    relcomponents_VaRs = components_VaRs / VaR

    return relcomponents_VaRs

'''#############################################################################################################################################################################'''
