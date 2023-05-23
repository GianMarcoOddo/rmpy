

from typing import Union, List
import numpy as np

'''##################################################################################### SINGLE AND PORT NONPARAMETRIC VAR ###############################################################################################'''

## 1. yf_npVaR_single #########################################################################################################################################################
############################################################################################################################################################################

def yf_npVaR_single(ticker : str,
                    position : Union[int,float],
                    start_date : str, 
                    end_date: str = "today",
                    freq: str = "daily",
                    alpha: float = 0.01,
                    kind: str = "abs",
                    display: bool = True) -> Union[float, None]:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function yf_npVaR_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE asset position using historical prices obtained from Yahoo Finance (yfinance).
    
    The arguments of this function are:

    "ticker": the asset symbol or ticker of the company for which we want to calculate the VaR.

    "position": the size of the position (value).

    "start_date": the starting date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object.

    "end_date": the ending date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object. By default, this is set to "today" 
    which will use the current date as the ending date.

    "freq": the frequency of the historical prices. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".

    "alpha": the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).

    "kind": the type of VaR calculation to use. This can be set to "abs" for absolute VaR or "rel" for relative VaR. By default, this is set to "abs".

    "display": a boolean value or string representing whether to display the calculated VaR. This can be set to True or False. By default, this is set to True.

    *******************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import yf_npVaR_single
    >>> ticker = ['AAPL'] 
    >>> position = [-1000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> kind = "abs"
    >>> display = True
    >>> yf_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha,kind = kind, display=display)
    >>> # OR
    >>> VaR = yf_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha,kind = kind, display=False)
    >>> print(VaR)

    This example calculates the daily non-parametric VaR for a short position of 1000 in Apple Inc, with a confidence level of 99%.

    *******************************************************************************************************************************************************************************************
    '''
    
    import numpy as np ; import pandas as pd ; import yfinance as yf
    from rmpy.inputs_handler import check_ticker_single, check_position_single , check_and_convert_dates , check_freq , check_kind , check_alpha , check_display

    # check ticker
    ticker = check_ticker_single(ticker, (str, list, np.ndarray, pd.Series), "ticker", single_value=True)
    # check position
    position = check_position_single(position) 
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check kind
    kind = check_kind(kind, allowed_values=["abs", "rel"])
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

    from rmpy.NpVaR import npVaR_single
    
    VaR = npVaR_single(returns, position, alpha = alpha, kind= kind)

    # Print settings
    freq_to_num = {"daily": 1,"weekly": 5,"monthly": 25}
    num = freq_to_num.get(freq)
    conf_level = int((1-alpha)*100)
    display = str(display).lower() in ['true', 't']
    if display:
            print(f"{kind.capitalize()} Nonparametric VaR({conf_level}%,{num}) for {ticker} =  {max(0,VaR):.3f} $")
    return VaR


## 2. yf_npVaR_port ##########################################################################################################################################################
############################################################################################################################################################################

def yf_npVaR_port(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  freq: str = "daily",
                  alpha: float = 0.01,
                  kind: str = "abs",
                  display: bool = True) -> Union[float, None]:
    '''
    *********************************************************************************************************************************************************************************************************

    DEFINITION: Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *******************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function 'yf_npVaR_port' is designed to calculate the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using historical prices obtained from Yahoo Finance (yfinance).
    
    The arguments of this function are:

    "tickers": A list of strings representing the tickers of the assets in the portfolio. #### note that all the TICKERS provided should be part of the portfolio whose VaR is being calculated ####

    "positions": A list of integers or floats representing the positions of the assets in the portfolio. The length of this list should be the same as the 'tickers' list.

    "start_date": A string representing the start date for the historical price data in the format 'YYYY-MM-DD'.

    "end_date": A string representing the end date for the historical price data in the format 'YYYY-MM-DD'. By default, it is set to "today".

    "freq": A string representing the frequency of the price data. By default, it is set to "daily".

    "alpha": A float representing the confidence level for the VaR calculation. By default, it is set to 0.01.

    "kind": A string representing the type of VaR calculation to perform. It can be either "abs" for absolute VaR or "rel" for relative VaR. By default, it is set to "abs".

    "display": A boolean indicating whether to print the VaR results or not. By default, it is set to True.

    *********************************************************************************************************************************************************************************************************
    
    Example
    -------
    >>> from rmpy.NpVaR import yf_npVaR_port
    >>> tickers = ['AAPL','MSFT'] 
    >>> positions = [-1000,5000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> kind = "abs"
    >>> display = True
    >>> yf_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha,kind = kind, display=display)
    >>> # OR
    >>> VaR = yf_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha,kind = kind, display=False)
    >>> print(VaR)

    This example calculates the daily non-parametric VaR for a Portfolio with short position of 1000 in Apple Inc and a long position of 5000 in Microsoft Corp.

    *****************************************************************************************************************************************************************************************************
    '''
    
    import yfinance as yf
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_freq , check_kind , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check freq
    freq = check_freq(freq, allowed_values=["daily", "weekly", "monthly"])
    # check alpha
    alpha = check_alpha(alpha)
    # check kind
    kind = check_kind(kind, allowed_values=["abs", "rel"])
    # check display
    check_display(display, allowed_values=[True, "True", "T", False, "False", "F"])

    # DOWNLOADING PROCESS  
    
    def get_data(tickers, start_date, end_date, freq, kind):
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
            raise Exception(f"{tickers_str} is/are incorrect or may have been delisted. Please check if this/these ticker/s is/are correct")
        return returns
    
    returns = get_data(tickers, start_date, end_date, freq, kind)

    from rmpy.NpVaR import npVaR_port
    
    VaR = npVaR_port(returns, positions, alpha = alpha,kind= kind)

    # Print settings
    freq_to_num = {"daily": 1,"weekly": 5,"monthly": 25}
    num = freq_to_num.get(freq)
    conf_level = int((1-alpha)*100)
    display = str(display).lower() in ['true', 't']
    if display:
            print(f"{kind.capitalize()} Nonparametric VaR({conf_level}%,{num}) for the Portfolio =  {VaR:.3f} $")
    return VaR

        
## 3. npVaR_single ############################################################################################################################################################
############################################################################################################################################################################

def npVaR_single(returns: np.ndarray,
                position : Union[int,float],
                 alpha: float = 0.01,
                 kind: str = "abs") -> Union[float, None]:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function npVaR_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE asset using historical returns data (you can use every type of assets e.g stock,options,bonds, ecc.)
    
    The arguments of this function are:

    "returns": a pandas Series or NumPy array containing the historical returns of the asset.

    "position": the size of the position in units of the asset. 

    "alpha": the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).

    "kind": the type of VaR calculation to use. This can be set to "abs" for absolute VaR or "rel" for relative VaR. By default, this is set to "abs".

    The function checks the "kind" argument to determine whether to calculate an absolute VaR or a relative VaR. If kind is set to "abs", the function calculates the absolute VaR, 
    which is the dollar amount of the potential loss. If kind is set to "rel", the function calculates the relative VaR, which is the percentage amount of the potential loss relative to the position size.
    If position is greater than 0, the function calculates the VaR for a long position. If position is less than 0, the function calculates the VaR for a short position. If position is 0, the function returns 
    a VaR of 0.
    The function first calculates the empirical quantile of the returns distribution at the specified confidence level (alpha). It then calculates the VaR based on the position size and the quantile value. 
    For absolute VaR, the VaR is the negative of the quantile value multiplied by the position size. For relative VaR, the VaR is the negative of the difference between the quantile value and the mean return, 
    multiplied by the position size.
    The function then returns the VaR value, but with a lower bound of 0. This ensures that the VaR is never negative

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import npVaR_single ; import numpy as np
    >>> returns = np.array([-0.02, 0.03, -0.01, 0.015, -0.002, 0.001, 0.008, 0.002, -0.006, 0.009]) # Replace this with actual returns
    >>> position = [-1000]
    >>> alpha = 0.01
    >>> kind = "abs"
    >>> VaR = npVaR_single(returns, position, alpha=alpha, kind = kind)
    >>> print(VaR)

    This example calculates the absolute non-parametric VaR consisting of a single short position of 1000 for the given returns.

    *****************************************************************************************************************************************************************************************************************
    '''
    import pandas as pd ; import numpy as np
    from rmpy.inputs_handler import validate_returns_single , check_position_single, check_alpha, check_kind

    # check returns
    returns = validate_returns_single(returns)
    # check position
    position = check_position_single(position) 
    # check alpha
    alpha = check_alpha(alpha)
    # check kind
    kind = check_kind(kind, allowed_values=["abs", "rel"])

    if kind == "abs":
        if position > 0: 
            quantile_alpha = np.quantile(returns, alpha)
            VaR = position * - quantile_alpha 
        if position < 0: 
            quantile_alpha = np.quantile(returns, 1- alpha)
            VaR = - position * quantile_alpha 
        if position == 0:
            VaR = 0
    if kind == "rel":
        mu = np.mean(returns)
        if position > 0: 
            quantile_alpha = np.quantile(returns, alpha)
            VaR = position * ( - quantile_alpha + mu )
        if position < 0: 
            quantile_alpha = np.quantile(returns, 1- alpha)
            VaR = - position * ( quantile_alpha - mu )
        if position == 0:
            VaR = 0

    return max(0,VaR)


## 4. npVaR_port ##############################################################################################################################################################
############################################################################################################################################################################

def npVaR_port(returns: np.ndarray,
            positions: List[Union[int, float]],
            alpha: float = 0.01,
            kind: str = "abs") -> Union[float, None]:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function npVaR_port calculates the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using historical returns data (you can use every type of assets e.g stock, options, bonds, ecc.) 
    
    The arguments of this function are:

    "returns": a pandas Series or NumPy array containing the historical returns of the portfolio.  #### note that all the RETURNS provided should be part of the portfolio whose VaR is being calculated ####

    "positions": the size of the positionS in units of the portfolio. This can be a single value or an array of values corresponding to each element in the returns argument.

    "alpha": the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).

    "kind": the type of VaR calculation to use. This can be set to "abs" for absolute VaR or "rel" for relative VaR. By default, this is set to "abs".

    Inside the function, the "kind" argument is checked to determine whether to calculate an absolute VaR or a relative VaR. The function then calculates the VaR based on the positions and returns of the portfolio.
    For absolute VaR, the function creates a copy of the returns DataFrame and iterates over each time period to calculate the portfolio value. It then calculates the VaR based on the quantile of the portfolio returns 
    distribution at the specified confidence level. For relative VaR, the function follows the same process as for absolute VaR, but subtracts the mean portfolio return from the quantile value before multiplying by the 
    position size.
    The function then returns the VaR value, but with a lower bound of 0. This ensures that the VaR is never negative.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import npVaR_port ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [-1000,2500,7000]
    >>> alpha = 0.01
    >>> kind = "abs"
    >>> VaR = npVaR_port(returns, position, alpha=alpha, kind = kind)
    >>> print(VaR)

    This example calculates the relative non-parametric VaR consisting of a portfolio with short positions of 1000 in the first asset, long positions of 2500 in the second asset, 
    and long positions of 7000 in the third asset.
    
    *************************************************************************************************************************************************************************************************************************
    '''
    import numpy as np
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_alpha, check_kind

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check alpha
    alpha = check_alpha(alpha)
    # check kind
    kind = check_kind(kind, allowed_values=["abs", "rel"])

        
    if kind == "abs":
        df_abs = np.copy(returns)
        # Every Position in the Positions array must be a float number  
        positions_float = []
        for num in positions:
            if type(num) == int:
                positions_float.append(float(num))
            else:
                positions_float.append(num)
        # Starting point --- this part calculates the evolution of the Porflolio over time
        new_positions = positions_float
        for i in range(df_abs.shape[0]):
            df_abs[i, :] *= new_positions 
            new_positions += df_abs[i, :]   
        port_rets = df_abs.sum(axis=1)
        # VaR
        VaR = np.quantile(port_rets,alpha) * -1
        return VaR
    
    if kind == "rel":
        df_rel = np.copy(returns)
        positions_float = []
        for num in positions:
            if type(num) == int:
                positions_float.append(float(num))
            else:
                positions_float.append(num)
        # Starting point --- this part calculates the evolution of the Porflolio over time
        new_positions = positions_float
        for i in range(df_rel.shape[0]):
            df_rel[i, :] *= new_positions
            new_positions += df_rel[i, :]
        port_rets = df_rel.sum(axis=1)
        # VaR
        mean_ret = np.mean(port_rets)
        VaR = (np.quantile(port_rets,alpha) - mean_ret) * -1
        
        return max(0,VaR)
    
'''########################################################## CONFIDENCE LEVEL SINGLE AND PORT NONPARAMETRIC VAR #############################################################'''

## 5. yf_conflevel_npVaR_single ###############################################################################################################################################
############################################################################################################################################################################

def yf_conflevel_npVaR_single(ticker : str,
                    position : Union[int,float],
                    start_date : str, 
                    end_date: str = "today",
                    freq: str = "daily",
                    alpha: float = 0.01,
                    confidence_level: float = 0.95,
                    display: bool = True) -> dict:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function yf_conflevel_npVaR_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE asset using Yahoo Finance data. It first downloads historical price data from Yahoo Finance, 
    calculates the returns of the asset, and then calculates the VaR a its confidence level (lower and upper bound)

    The arguments of this function are:

    "ticker": the symbol of the asset to calculate VaR for.

    "position": the size of the position in units of the asset.

    "start_date": the start date of the historical data to download. This should be a string in the format "YYYY-MM-DD".

    "end_date": the end date of the historical data to download. This should be a string in the format "YYYY-MM-DD". By default, this is set to "today".

    "freq": the frequency of the data to download. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".

    "alpha": the level of significance for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
    
    "confidence_level": the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability that the true VaR is within the calculated VaR interval. 
    By default, this is set to 0.95 (95%).

    "display": a boolean value indicating whether to display the VaR calculation result. By default, this is set to True.

    The function then, is display is true, prints the VaR calculation with its confidence interval.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import yf_conflevel_npVaR_single
    >>> ticker = ['AAPL']
    >>> position = [2500]
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> confidence_level = 0.95
    >>> display = True
    >>> yf_conflevel_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha ,confidence_level = confidence_level, display=display)
    >>> # OR
    >>> VaR = yf_conflevel_npVaR_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha ,confidence_level = confidence_level, display=False)
    >>> print(VaR)

    This example calculates the absolute non-parametric VaR with its lower and upper bound (95% confidence) for a position of 2500 in "AAPL".

    *********************************************************************************************************************************************************************************************************
    '''
    import numpy as np ; import pandas as pd ; import yfinance as yf
    from rmpy.inputs_handler import check_ticker_single, check_position_single , check_and_convert_dates , check_freq , check_confidence_level , check_alpha , check_display

    # check ticker
    ticker = check_ticker_single(ticker, (str, list, np.ndarray, pd.Series), "ticker", single_value=True)
    # check position
    position = check_position_single(position) 
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
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
            raise Exception(f"Insufficient data to calculate the Nonparametric VaR for {ticker}")
        return returns

    # Call function with your parameters
    returns = get_data(ticker, start_date, end_date, freq)

    from rmpy.NpVaR import conflevel_npVaR_single

    VaR = conflevel_npVaR_single(returns,position,confidence_level,alpha)
    # Print settings
    display = str(display).lower() in ['true', 't']
    if display:
            rounded_VaR = {key: round(value, 4) for key, value in VaR.items()}
            print(rounded_VaR)
    return VaR

## 6. yf_conflevel_npVaR_port ####################################################################################################################################################
############################################################################################################################################################################

def yf_conflevel_npVaR_port(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  freq: str = "daily",
                  alpha: float = 0.01,
                  confidence_level : float = 0.95,
                  display: bool = True) -> dict:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function yf_conflevel_npVaR_port calculates the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using Yahoo Finance data. It first downloads historical price data from Yahoo Finance
    for each asset in the portfolio, calculates the returns of each asset, and then  calculates the portfolio VaR, the lower and the upper bound at a specified confidence level and alpha level.

    The arguments of this function are:

    "tickers": a list of symbols of the assets in the portfolio to calculate VaR for.

    "positions": a list or array containing the sizes of the positions for each asset in the portfolio. These sizes are in units of the assets and can be positive or negative.

    "start_date": the start date of the historical data to download. This should be a string in the format "YYYY-MM-DD".

    "end_date": the end date of the historical data to download. This should be a string in the format "YYYY-MM-DD". By default, this is set to "today".

    "freq": the frequency of the data to download. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".

    "alpha": the level of significance for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
   
    "confidence_level": the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability that the true VaR is within the calculated VaR interval. By default, this is set to 0.95.
    
    "display": a boolean value indicating whether to display the VaR calculation result. By default, this is set to True.

    The function first imports the conflevel_npVaR_port function from the rmpy.NpVaR module, which calculates the portfolio VaR and its cofidence interval using the non-parametric method at a specified confidence level and alpha level.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import yf_conflevel_npVaR_port
    >>> tickers = ['AAPL','MSFT'] 
    >>> positions = [-1000,5000] 
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> confidence_level = 0.95
    >>> display = True
    >>> yf_conflevel_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha,confidence_level=confidence_level, display=display)
    >>> # OR
    >>> VaR = yf_conflevel_npVaR_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha,confidence_level=confidence_level, display=False)
    >>> print(VaR)

    This example calculates the daily non-parametric VaR with its lower and upper bound (95% confidence) for a Portfolio with short position of 1000 in Apple Inc and 5000 in Microsoft Corp.

    ***********************************************************************************************************************************************************************************************************
    '''
    import yfinance as yf 
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_freq , check_confidence_level , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
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
        # If the provided frequency is not supported, raise an exception
        if freq not in freq_to_interval:
            raise ValueError(f"Unsupported frequency '{freq}'")
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

    from rmpy.NpVaR import conflevel_npVaR_port
    
    VaR = conflevel_npVaR_port(returns, positions, confidence_level, alpha)
    # Print settings
    display = str(display).lower() in ['true', 't']
    if display:
            rounded_VaR = {key: round(value, 4) for key, value in VaR.items()}
            print(rounded_VaR)
    return VaR

## 7. conflevel_npVaR_single ##################################################################################################################################################
############################################################################################################################################################################

def conflevel_npVaR_single(returns: np.ndarray,
                           position: Union[int, float],
                           confidence_level: float = 0.95,
                           alpha: float = 0.01)-> dict:

    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function conflevel_npVaR_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE asset using historical returns data (you can use every type of assets e.g stock,options,bonds, ecc.) with 
    a specific confidence interval and alpha value. It also caluculates the lower and upper bound for the non-parametric VaR.

    The arguments of this function are:

    "returns": A NumPy array or pd.Series of historical returns of a single asset.

    "position": The size of the position in the asset. If position is greater than 0, the function assumes the position is long, and if position is less than or equal to 0, the function assumes the position is short.
    
    "confidence_level": The confidence level for the VaR calculation (default value is 0.95).
    
    "alpha": The level of significance for the VaR calculation (default value is 0.01).

    The function first checks whether the position is long or short, and then calculates the VaR using the appropriate formula for the chosen confidence level and alpha.
    If the position is long, the function uses the upper quantile of the returns distribution to calculate VaR, while if the position is short, the function uses the lower quantile.
    The output of the function is a dictionary containing three keys: npvar, lower_bound, and upper_bound. npvar is the Non-Parametric Value at Risk, lower_bound is the lower bound of the confidence interval, 
    and upper_bound is the upper bound of the confidence interval. The values of lower_bound and upper_bound are calculated using the standard error of the quantile, the inverse cumulative distribution function of 
    either the normal or the t-student distribution in automatic, and the position size.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import conflevel_npVaR_single ; import numpy as np
    >>> returns = np.random.uniform(-0.03, 0.03, size=1000) # Replace this with actual returns
    >>> position = [-1000]
    >>> confidence_level = 0.95
    >>> alpha = 0.01
    >>> VaR = conflevel_npVaR_single(returns, position,confidence_level= confidence_level, alpha=alpha)
    >>> print(VaR)

    This example calculates the non-parametric VaR and its lower and upper bound (95% confidence interval) consisting of a position of -1000 in the given asset.

    *****************************************************************************************************************************************************************************************************************
    '''
    import pandas as pd ; import numpy as np ;  from scipy import stats ; import warnings
    from rmpy.inputs_handler import validate_returns_single , check_position_single, check_confidence_level, check_alpha

    # check returns
    returns = validate_returns_single(returns)
    # check position
    position = check_position_single(position) 
    # check confidence_level 
    confidence_level = check_confidence_level(confidence_level)
    # check alpha
    alpha = check_alpha(alpha)

    # Testing for Normality

    # Suppress the warning message
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Perform the Shapiro-Wilk test
        W, p = stats.shapiro(returns)

    # LONG #################################################################

    if position > 0:

        if p > 0.05:
            mu = np.mean(returns)
            sigma = np.std(returns)
            from scipy.stats import norm
            quantile = np.quantile(returns,alpha) 
            se_q = (1 / norm.pdf(quantile, loc=mu, scale=sigma) ) * np.sqrt ( ( ( 1-alpha) * alpha ) /   returns.shape[0]  )
            alpha_conf_level = norm.ppf((1-confidence_level)/2) * -1
            quantile = quantile * -1
            var = quantile * position
            if var <= 0:
                var = 0 
            lower_bound = quantile - se_q * alpha_conf_level
            upper_bound = quantile + se_q * alpha_conf_level
            VaR = {"npvar": var, "lower_bound":position * lower_bound,"upper_bound": position * upper_bound}
            return VaR
        
        if p <= 0.05:
            from scipy.stats import t
            from scipy.stats import norm
            df, loc, scale = t.fit(returns) 
            quantile = np.quantile(returns,alpha)
            se_q = (1 / t.pdf(quantile, df, loc, scale)) * np.sqrt ( ( ( 1-alpha) * alpha ) /   returns.shape[0]  )
            alpha_conf_level = norm.ppf((1-confidence_level)/2) * -1
            quantile = quantile * -1
            var = quantile * position
            if var <= 0:
                var = 0 
            lower_bound = quantile - se_q * alpha_conf_level
            upper_bound = quantile + se_q * alpha_conf_level
            VaR = {"npvar": var, "lower_bound": position * lower_bound,"upper_bound": position * upper_bound}
            return VaR

            
    # SHORT #################################################################
            
    if position <= 0:
        if p > 0.05:
            mu = np.mean(returns)
            sigma = np.std(returns)
            from scipy.stats import norm
            quantile = np.quantile(returns,1-alpha)
            se_q = (1 / norm.pdf(quantile, loc=mu, scale=sigma) ) * np.sqrt ( ( ( 1-alpha) * alpha ) /   returns.shape[0])
            alpha_conf_level = norm.ppf((1-confidence_level)/2) * -1
            var = quantile * - position
            if var < 0:
                var = 0
            lower_bound = quantile - se_q * alpha_conf_level
            upper_bound = quantile + se_q * alpha_conf_level
            VaR = {"npvar": var, "lower_bound":  - position * lower_bound,"upper_bound": - position * upper_bound}
            return VaR

        if p <= 0.05:
            from scipy.stats import t
            from scipy.stats import norm
            df, loc, scale = t.fit(returns)
            quantile = np.quantile(returns,1-alpha)
            se_q = (1 / t.pdf(quantile, df, loc, scale)) * np.sqrt ( ( ( 1-alpha) * alpha ) /   returns.shape[0]  )
            alpha_conf_level = norm.ppf((1-confidence_level)/2) * -1
            var = quantile * - position
            if var < 0:
                var = 0
            lower_bound = quantile - se_q * alpha_conf_level
            upper_bound = quantile + se_q * alpha_conf_level
            VaR = {"npvar": var, "lower_bound":  - position * lower_bound,"upper_bound": - position * upper_bound}
            return VaR
    

## 8. conflevel_npVaR_port ####################################################################################################################################################
############################################################################################################################################################################

def conflevel_npVaR_port(returns : np.ndarray,
                        positions: List[Union[int, float]],
                        confidence_level: float = 0.95,
                        alpha: float = 0.01)-> dict:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function conflevel_npVaR_port calculates the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using historical returns data (you can use every type of assets e.g stock,options,bonds, ecc.) with 
    a specific confidence interval and alpha value. It also calculates the lower and the upper bound of this estimation. 

    The arguments of this function are:

    "returns": A NumPy array of historical returns for the portfolio.
    
    "positions": A NumPy array of positions for each asset in the portfolio.
    
    "confidence_level": The confidence level for the VaR calculation (default value is 0.95).
    
    "alpha": The level of significance for the VaR calculation (default value is 0.01).

    The function first calculates the portfolio returns by multiplying the historical returns by the positions for each asset,
    and then sums the products to obtain the portfolio returns.
    The output of the function is a dictionary containing three keys: npvar, lower_bound, and upper_bound. npvar is the Non-Parametric Value at Risk, lower_bound is the lower bound of the confidence interval, 
    and upper_bound is the upper bound of the confidence interval.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import conflevel_npVaR_port ; import numpy as np
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [-1000,2000,-3000]
    >>> confidence_level = 0.95
    >>> alpha = 0.01
    >>> VaR = conflevel_npVaR_port(returns, positions ,confidence_level= confidence_level, alpha=alpha)
    >>> print(VaR)

    This example calculates the non-parametric VaR consisting and its lower and upper bound (95% CI) of portfolio consisting of a postion of -1000 in the first asset,
    2000 in the second asset and - 3000 in the third one.

    ********************************************************************************************************************************************************************************************************************************************
    '''
    import numpy as np ; from scipy import stats ; import warnings
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_confidence_level, check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check position
    position = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check confidence_level 
    confidence_level = check_confidence_level(confidence_level)
    # check alpha
    alpha = check_alpha(alpha)
    
    # Every Position in the Positions array must be a float number
    positions_float = []
    for num in positions:
        if type(num) == int:
            positions_float.append(float(num))
        else:
            positions_float.append(num)
    df_temp = np.copy(returns)

    # Starting point --- this part calculates the evolution of the Porflolio over time

    new_positions = positions_float
    for i in range(df_temp.shape[0]):
        df_temp[i, :] *= new_positions
        new_positions += df_temp[i, :]

    port_rets = df_temp.sum(axis=1)

    # Testing for Normality

    # Suppress the warning message
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Perform the Shapiro-Wilk test
        W, p = stats.shapiro(returns)

    if p > 0.05:
        mu = np.mean(port_rets)
        sigma = np.std(port_rets)
        from scipy.stats import norm
        quantile = np.quantile(port_rets,alpha)
        se_q = (1 / norm.pdf(quantile, loc=mu, scale=sigma) ) * np.sqrt ( ( ( 1-alpha) * alpha ) /   port_rets.shape[0])
        alpha_conf_level = norm.ppf((1-confidence_level)/2) * -1
        var = - quantile 
        if var < 0:
            var = 0
        lower_bound = - quantile - se_q * alpha_conf_level
        upper_bound = - quantile + se_q * alpha_conf_level
        VaR = {"npvar": var, "lower_bound": lower_bound,"upper_bound": upper_bound}
        return VaR

    if p <= 0.05:
        from scipy.stats import t
        from scipy.stats import norm
        df, loc, scale = t.fit(port_rets)
        quantile = np.quantile(port_rets,alpha)
        se_q = (1 / t.pdf(quantile, df, loc, scale)) * np.sqrt ( ( ( 1-alpha) * alpha ) /   port_rets.shape[0]  )
        alpha_conf_level = norm.ppf((1-confidence_level)/2) * -1
        var = - quantile 
        if var < 0:
            var = 0
        lower_bound = - quantile - se_q * alpha_conf_level
        upper_bound = - quantile + se_q * alpha_conf_level
        VaR = {"npvar": var, "lower_bound": lower_bound,"upper_bound": upper_bound}
        return VaR
        

'''############################################################### # SUMMARY - SINGLE AND PORTFOLIO  ######################################################################'''

## 9. yf_npVaR_summary_single #################################################################################################################################################
############################################################################################################################################################################  

def yf_npVaR_summary_single(ticker: str, 
                            position: Union[int, float],
                            start_date: str,
                            end_date: str = "today",
                            freq: str = "daily",
                            alpha: float = 0.01,
                            display: bool = True)-> dict:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function yf_npVaR_summary_single calculates the quantile non-parametric Value at Risk (VaR) for a SINGLE stock position using historical prices obtained from Yahoo Finance (yfinance).
    This function is usefull to visualize some risk metrics of the single asset npVaR.

    The arguments of this function are:

    "ticker": the stock symbol or ticker of the company for which we want to calculate the VaR. This can be a string (for a single ticker) or a list or array of tickers (for multiple tickers).
    
    "position": the size of the position in shares or currency. This can be a single value or an array of values corresponding to each ticker in the ticker argument.
    
    "start_date": the starting date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object.
    
    "end_date": the ending date for which we want to obtain historical prices. This can be a string in the format "YYYY-MM-DD" or a datetime object. By default, this is set to "today" 
     which will use the current date as the ending date.
    
    "freq": the frequency of the historical prices. This can be set to "daily", "weekly", or "monthly". By default, this is set to "daily".
    
    "alpha": the confidence level for the VaR calculation. This is a value between 0 and 1, representing the probability of loss exceeding the VaR. By default, this is set to 0.01 (1%).
    
    "display": a boolean value or string representing whether to display the calculated VaR. This can be set to True, "True", "T", False, "False", or "F". By default, this is set to True.

    The function prints the summary of the VaR calculation if display is True and returns the summary as a dictionary. 
    The summary contains the following key-value pairs:

    npVaR: The non-parametric Value at Risk (VaR)

    max_loss: The maximum loss

    max_excess_loss: The maximum excess loss

    max_excess_loss_over_npVaR: The maximum excess loss over the non-parametric VaR

    expected_shortfall: The expected shortfall

    expected_shortfall_over_npVaR: The expected shortfall over the non-parametric VaR

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import yf_npVaR_summary_single
    >>> ticker = 'AAPL'
    >>> position = -1000
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_npVaR_summary_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha, display=display)
    >>> # OR 
    >>> VaR = yf_npVaR_summary_single(ticker, position, start_date, end_date, freq=freq, alpha=alpha, display=False)

    This example calculates several metrics for the non-parametric VaR summary for a short position of 1000 shares in Apple Inc with a confidence level of 99% using historical daily prices
    from January 1, 2020, to December 31, 2021, obtained from Yahoo Finance.

    *******************************************************************************************************************************************************************************************
    '''
    import numpy as np ; import pandas as pd ; import yfinance as yf 
    from rmpy.inputs_handler import check_ticker_single, check_position_single , check_and_convert_dates , check_freq  , check_alpha , check_display

    # check ticker
    ticker = check_ticker_single(ticker, (str, list, np.ndarray, pd.Series), "ticker", single_value=True)
    # check position
    position = check_position_single(position) 
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
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
        if data.empty:
            raise Exception(f"The provided ticker {ticker} may have been delisted or it might be wrong")

        # Calculate returns and drop missing values
        returns = data["Adj Close"].pct_change().dropna()

        # Raise an exception if less than two returns have been downloaded
        if returns.shape[0] < 2:
            raise Exception(f"Insufficient data to calculate the Nonparametric VaR for {ticker}")
        return returns

    # Call function with your parameters
    returns = get_data(ticker, start_date, end_date, freq)

    from rmpy.NpVaR import npVaR_summary_single

    summary = npVaR_summary_single(returns, position, alpha =alpha)
    # Print settings
    display = str(display).lower() in ['true', 't']
    if display:
            rounded_summary = {key: round(value, 4) for key, value in summary.items()}
            print(rounded_summary)
    return summary


## 10. yf_npVaR_summary_port  ##################################################################################################################################################
############################################################################################################################################################################  

def yf_npVaR_summary_port(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  freq: str = "daily",
                  alpha: float = 0.01,
                  display: bool = True) -> dict:
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function 'yf_npVaR_summary_port' is designed to calculate the quantile non-parametric Value at Risk (VaR) for a PORTFOLIO of assets using historical prices obtained from Yahoo Finance (yfinance).
    This function is usefull to visualize several metrics of the Portfolio npVaR.

    The argument of this function are:

    "tickers": A list of strings representing the tickers of the assets in the portfolio. #### note that all the TICKERS provided should be part of the portfolio whose VaR is being calculated ####
    
    "positions": A list of integers or floats representing the positions of the assets in the portfolio. The length of this list should be the same as the 'tickers' list.
    
    "start_date": A string representing the start date for the historical price data in the format 'YYYY-MM-DD'.
    
    "end_date": A string representing the end date for the historical price data in the format 'YYYY-MM-DD'. By default, it is set to "today".
    
    "freq": A string representing the frequency of the price data. By default, it is set to "daily".
    
    "alpha": A float representing the confidence level for the VaR calculation. By default, it is set to 0.01.
    
    "kind": A string representing the type of VaR calculation to perform. It can be either "abs" for absolute VaR or "rel" for relative VaR. By default, it is set to "abs".
   
    "display": A boolean indicating whether to print the VaR results or not. By default, it is set to True.

    The function prints the summary of the VaR calculation if display is True and returns the summary as a dictionary. 
    The summary contains the following key-value pairs:

    npVaR: The non-parametric Value at Risk (VaR)
    max_loss: The maximum loss
    max_excess_loss: The maximum excess loss
    max_excess_loss_over_npVaR: The maximum excess loss over the non-parametric VaR
    expected_shortfall: The expected shortfall
    expected_shortfall_over_npVaR: The expected shortfall over the non-parametric VaR

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import yf_npVaR_summary_port
    >>> tickers = ['AAPL', 'MSFT']
    >>> positions = [-1000, 5000]
    >>> start_date = '2020-01-01'
    >>> end_date = '2021-12-31'
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_npVaR_summary_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, display=display)
    >>> # OR 
    >>> VaR = yf_npVaR_summary_port(tickers, positions, start_date, end_date, freq=freq, alpha=alpha, display=False)

    This example calculates several metrics for the non-parametric VaR for a portoflio consisting of a short position of 1000 in Apple Inc and a long position of
    5000 in Miscrosoft using historical daily prices from January 1, 2020, to December 31, 2021, obtained from Yahoo Finance.

    *********************************************************************************************************************************************************************************************************
    '''
    import yfinance as yf
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates , check_freq , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
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
        # If the provided frequency is not supported, raise an exception
        if freq not in freq_to_interval:
            raise ValueError(f"Unsupported frequency '{freq}'")
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

    from rmpy.NpVaR import npVaR_summary_port

    summary = npVaR_summary_port(returns, positions, alpha)
    # Print settings
    display = str(display).lower() in ['true', 't']
    if display:
            rounded_summary = {key: round(value, 4) for key, value in summary.items()}
            print(rounded_summary)
    return summary


## 11. npVaR_summary_single ####################################################################################################################################################
############################################################################################################################################################################  

def npVaR_summary_single(returns: np.ndarray,
                         position: Union[int, float],
                        alpha: float = 0.01)-> dict:

    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    This function called npVaR_summary_single that calculates the quantile Non-Parametric Value at Risk (VaR) and several metrics for a SINGLE asset and returns a summary of the VaR calculation.

    The argument of this function are:

    "returns": A numpy array or pandas series of historical returns for the asset.

    "position": The size of the position in the asset.

    "alpha": The level of significance for the VaR calculation. Default value is 0.01.

    The function first calculates the VaR for the single asset using the Non-Parametric method. If the position is positive, the VaR is calculated as the negative of the product of the position and the alpha-quantile of the returns. If the position is negative or zero, the VaR is calculated as the product of the position and the 1-alpha-quantile of the returns.
    Then, the function calculates additional risk measures using the VaR. The additional risk measures include: maximum loss, maximum excess loss, maximum excess loss over VaR, expected shortfall, and expected shortfall over VaR.
    
    The function returns the summary as a dictionary containing the following keys and corresponding values:

    npVaR: The Non-Parametric Value at Risk (npVaR) for the asset.

    max_loss: The maximum loss that could be incurred in the worst-case scenario.

    max_excess_loss: The maximum loss that exceeds the VaR.

    max_excess_loss_over_npVaR: The ratio of maximum excess loss over the VaR.

    expected_shortfall: The expected loss given that the return falls below the VaR.

    expected_shortfall_over_npVaR: The ratio of expected shortfall over the VaR.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> import numpy as np
    >>> from rmpy.NpVaR import npVaR_summary_single
    >>> returns = np.array([-0.02, 0.03, -0.01, 0.015, -0.002, 0.001, 0.008, 0.002, -0.006, 0.009]) # Replace this with actual returns
    >>> position = 10000
    >>> alpha = 0.01
    >>> summary = npVaR_summary_single(returns, position, alpha = 0.01)
    >>> print(summary)

    This example calculates all the metrics for the non-parametric VaR for a single asset with a given set of historical returns and position.

    *********************************************************************************************************************************************************************************************************
    '''
    import pandas as pd ; import numpy as np
    from rmpy.inputs_handler import validate_returns_single , check_position_single, check_alpha

    # check returns
    returns = validate_returns_single(returns)
    # check position
    position = check_position_single(position) 
    # check alpha
    alpha = check_alpha(alpha)
    
    #################################################################################################################
    
    # VaR

    if position > 0:
        quantile_alpha = np.quantile(returns,alpha)
        npvar = position * quantile_alpha
        VaR = - npvar
        Max_loss = np.min(returns) * - position
        Max_Excess_loss = Max_loss - VaR
        Max_loss_over_VaR = Max_Excess_loss/VaR
        Expected_shortfall = np.mean(returns[returns<quantile_alpha]) * - position
        Expected_shortfall_over_VaR = Expected_shortfall/VaR
        summary = {"npVaR": round(VaR,4), "max_loss": round(Max_loss,4),"max_excess_loss": round(Max_Excess_loss,4),
                "max_excess_loss_over_npVaR":round(Max_loss_over_VaR,4),"expected_shortfall":round(Expected_shortfall,4),
                "expected_shortfall_over_npVaR":round(Expected_shortfall_over_VaR,4)}
        return summary

    if position <= 0:
        quantile_alpha = np.quantile(returns,1-alpha)
        npvar = -1 * position * quantile_alpha
        VaR =  npvar
        Max_loss = np.max(returns) * - position
        Max_Excess_loss = Max_loss - VaR
        Max_loss_over_VaR = Max_Excess_loss/VaR
        Expected_shortfall = np.mean(returns[returns>quantile_alpha]) * - position
        Expected_shortfall_over_VaR = Expected_shortfall/VaR 
        summary = {"npVaR": round(VaR,4), "max_loss": round(Max_loss,4),"max_excess_loss": round(Max_Excess_loss,4),
                "max_excess_loss_over_npVaR":round(Max_loss_over_VaR,4),"expected_shortfall":round(Expected_shortfall,4),
                "expected_shortfall_over_npVaR":round(Expected_shortfall_over_VaR,4)}
        return summary
    

## 12. npVaR_summary_port ######################################################################################################################################################
############################################################################################################################################################################  
    
def npVaR_summary_port(returns: np.ndarray,
                       positions: List[Union[int, float]],
                       alpha: float = 0.01)-> dict:
    
    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    The function npVaR_summary_port calculates the quantile non-parametric Value at Risk (VaR) and several metrics for a npVaR of a portfolio.
    
    The argument of this function are:

    "returns": A numpy array or pandas series of historical returns for the portfolio.

    "positions": A pandas DataFrame or numpy array containing the value held for each asset in the portfolio;

    "alpha": The confidence level for the VaR calculation (default value is 0.01).

    The function first calculates the evolution of the portfolio over time based on the provided positions. It then calculates the VaR, maximum loss, maximum excess loss, maximum excess loss over VaR, 
    expected shortfall, and expected shortfall over VaR, all of which are returned in a dictionary. The VaR and expected shortfall are calculated using a non-parametric method based on the historical returns 
    of the portfolio.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> import numpy as np
    >>> from rmpy import npVaR_summary_port
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 3)) # Replace this with actual returns
    >>> positions = [1000, 500, -1500]
    >>> alpha = 0.01
    >>> summary = npVaR_summary_port(returns, positions, alpha=alpha)
    >>> print(summary)

    This example shows how to use the npVaR_summary_port function with a 10x3 numpy array of historical returns and a list of positions for 3 assets in the portfolio. 
    The function calculates the non-parametric Value at Risk (VaR) and other risk measures, returning them in a dictionary.

    *************************************************************************************************************************************************************************************************************************
    '''

    import numpy as np
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_alpha

    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check alpha
    alpha = check_alpha(alpha)
        
    df_abs = np.copy(returns)
    # Every Position in the Positions array must be a float number  
    positions_float = []
    for num in positions:
        if type(num) == int:
            positions_float.append(float(num))
        else:
            positions_float.append(num)
    # Starting point --- this part calculates the evolution of the Porflolio over time
    new_positions = positions_float
    for i in range(df_abs.shape[0]):
        df_abs[i, :] *= new_positions 
        new_positions += df_abs[i, :]   
    port_rets = df_abs.sum(axis=1)

    # VaR
    VaR = np.quantile(port_rets,alpha) * -1
    Max_loss = np.min(port_rets) * -1 
    Max_Excess_loss = Max_loss - VaR
    Max_loss_over_VaR = Max_Excess_loss/VaR
    Expected_shortfall = np.mean(port_rets[port_rets<np.quantile(port_rets,alpha)]) * -1
    Expected_shortfall_over_VaR = Expected_shortfall/VaR
    summary = {"npVaR": round(VaR,4), "max_loss": round(Max_loss,4),"max_excess_loss": round(Max_Excess_loss,4),
            "max_excess_loss_over_npVaR":round(Max_loss_over_VaR,4),"expected_shortfall":round(Expected_shortfall,4),
            "expected_shortfall_over_npVaR":round(Expected_shortfall_over_VaR,4)}
    return summary

'''################################################################## MARGINALS NON PARAMETRIC VAR  ######################################################################'''
    

## 13. yf_marg_NpVaRs_scale_factor ################################################################################################################################################
############################################################################################################################################################################

def yf_marg_NpVaRs_scale_factor(tickers: List[str], 
                  positions: List[Union[int, float]],
                  start_date: str,
                  end_date: str = "today",
                  scale_factor: float = 0.1,
                  freq: str = "daily",
                  alpha: float = 0.01,
                  display: bool = True)-> list:

    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    The function yf_marg_NpVaRs_scale_factor calculates the Marginal Non-Parametric VaRs and their changes for a PORTFOLIO of assets over a given period of time, using a specified scale factor,using historical 
    prices obtained from Yahoo Finance (yfinance). New position = old one + (old one * scale_factor) in the same direction (- if short and + if long)

    The argument of this function are:

    "tickers": a list of tickers for the assets in the portfolio.

    "positions": a list of positions for each asset in the portfolio.

    "start_date": the start date of the analysis period.

    "end_date": the end date of the analysis period (optional, default is "today").

    "scale_factor": the scale factor used to calculate the Marginal Non-Parametric VaRs (optional, default is 0.1).

    "freq": the frequency of the data (optional, default is "daily").

    "alpha": the confidence level used in calculating the VaRs (optional, default is 0.01).

    "display": a boolean indicating whether to print the resulting VaR changes (optional, default is True).

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> from rmpy.NpVaR import yf_marg_NpVaRs_scale_factor
    >>> tickers = ["AAPL", "MSFT", "GOOG"]
    >>> positions = [1000, 1500, 2000]
    >>> start_date = "2021-01-01"
    >>> end_date = "2021-12-31"
    >>> scale_factor = 0.1
    >>> freq = "daily"
    >>> alpha = 0.01
    >>> display = True
    >>> yf_marg_NpVaRs_scale_factor(tickers, positions, start_date, end_date, scale_factor=scale_factor, freq=freq, alpha=alpha, display=display)
    >>> #OR 
    >>> NpVaR_changes = yf_marg_NpVaRs_scale_factor(tickers, positions, start_date, end_date, scale_factor=scale_factor, freq=freq, alpha=alpha, display=False)
    >>> print(NpVaR_changes)  

    This example demonstrates how to use the yf_marg_NpVaRs_scale_factor function with a list of three tickers,
    a list of positions for each asset, a start date, and an end date for the analysis period. 
    The function calculates the Marginal Non-Parametric VaRs and their changes for the given portfolio using historical prices from Yahoo Finance.

    *************************************************************************************************************************************************************************************************************************
    '''
    import yfinance as yf 
   
    from rmpy.inputs_handler import check_tickers_port, check_positions_port , check_and_convert_dates ,check_scale_factor, check_freq , check_alpha , check_display

    # check tickers
    tickers = check_tickers_port(tickers)
    # chech positions
    positions = check_positions_port(positions)
    # check tickers and positions
    if len(tickers) != len(positions):
        raise ValueError("The number of tickers != the number of positions")
    # check dates
    start_date, end_date = check_and_convert_dates(start_date, end_date)
    # check scale_factor
    scale_factor = check_scale_factor(scale_factor)
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
        # If the provided frequency is not supported, raise an exception
        if freq not in freq_to_interval:
            raise ValueError(f"Unsupported frequency '{freq}'")
        # Perform the download once using the relevant interval
        interval = freq_to_interval[freq]
        data = yf.download(tickers, start=start_date, end=end_date, interval=interval, progress=False, show_errors=False)
        # Check if the data is empty
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
    
    from rmpy.NpVaR import marg_NpVaRs_scale_factor
    NpVaR_changes =  marg_NpVaRs_scale_factor(returns, positions, scale_factor = scale_factor, alpha = alpha)
    if display:
            rounded_NpVaR_changes = [round(value, 4) for value in NpVaR_changes]
            var_dict = dict(zip(tickers, rounded_NpVaR_changes))
            print(var_dict)
    return NpVaR_changes
 

## 14. marg_NpVaRs_scale_factor ################################################################################################################################################
############################################################################################################################################################################

def marg_NpVaRs_scale_factor(returns: np.ndarray,
                            positions: List[Union[int,float]],
                            scale_factor: float = 0.1,
                            alpha: float = 0.01) -> list:

    '''
    ******************************************************************************************************************************************************************************************************************
    
    DEFINITION
    -------
    Non-parametric VaR (Value at Risk) is a risk management tool used to estimate the potential losses in the value of an asset or portfolio. Unlike parametric VaR, 
    which is calculated based on assumed probability distributions, non-parametric VaR is based on historical price movements and does not rely on any distribution assumptions.

    Non-parametric VaR is typically used when there is not enough data to estimate the underlying probability distribution or when the asset returns do not follow a known distribution. 
    It is also useful when dealing with extreme market events where the assumption of a normal distribution may not hold.

    *****************************************************************************************************************************************************************************************************************
    
    Args
    -------
    The function marg_NpVaRs_scale_factor calculates the marginal changes in the non-parametric Value at Risk (NpVaR) of a PORTFOLIO based on a specified scale factor. 
    New position = old one + (old one * scale_factor) in the same direction (- if short and + if long)

    The argument of this function are:

    "returns": A numpy array or pandas series of historical returns for the portfolio.

    "positions": A list or array containing the current positions of the assets in the portfolio.

    "scale_factor": A float representing the scale factor to be used for calculating the marginal changes in NpVaR.

    "alpha": A float representing the significance level for calculating the NpVaR.

    *********************************************************************************************************************************************************************************************************

    Example
    -------
    >>> import numpy as np
    >>> from rmpy.NpVaR import marg_NpVaRs_scale_factor
    >>> returns = np.random.uniform(-0.05, 0.05, size=(10, 5)) # Replace this with actual returns
    >>> positions = [1000, -5000, 2000, -1000, 1500]
    >>> scale_factor = 0.2
    >>> alpha = 0.01
    >>> NpVaR_changes = marg_NpVaRs_scale_factor(returns, positions, alpha=alpha)
    >>> print(NpVaR_changes)

    In this example, we have calculated the marginal changes in the non-parametric value at risk (NpVaR) for a portfolio consisting of 5 assets. 
    The function takes in historical returns data for the assets, the current positions of the assets, a scale factor, and a significance level (alpha). 
    The result is a list of the marginal changes in NpVaR for each asset, which can be useful for understanding the risk contributions of individual 
    positions within the portfolio and designing risk management strategies to mitigate potential losses.

    **************************************************************************************************************************************************************************************************************************
    '''
    import numpy as np
    from rmpy.inputs_handler import validate_returns_port , check_positions_port, check_scale_factor, check_alpha
 
    # check returns
    returns = validate_returns_port(returns)
    # check positions
    positions = check_positions_port(positions)
    # THE NUMBER OF POSITIONS SHOULD MATCH THE NUMBER OF ROWS IN THE ARRAY OF RETURNS
    if len(positions) != returns.shape[1]:
        raise ValueError("The number of assets != the number of positions")
    # check scale_factor
    scale_factor = check_scale_factor(scale_factor)
    # check alpha
    alpha = check_alpha(alpha)

    positions_float = []
    for num in positions:
        if type(num) == int:
            positions_float.append(float(num))
        else:
            positions_float.append(num)

    ####### VaR Calculator ####### - abs 

    def Var_calculator(returns, positions_float, alpha):
        df_abs = np.copy(returns) ;  new_positions = positions_float
        for i in range(df_abs.shape[0]):
            df_abs[i, :] *= new_positions
            new_positions += df_abs[i, :]
        port_rets = df_abs.sum(axis=1)

        VaR = np.quantile(port_rets,alpha) * -1
        if VaR <= 0:
            VaR = 0
        return VaR

    ###########################################################################################

    # The Nonparametric VaR is affected by changes when each position is increased by an amount proportional to 1/scale_factor, in accordance with its corresponding sign

    NpVaR_changes = []
    for index in range(len(positions_float)):
        prev_var = Var_calculator(returns, positions_float, alpha)
        if positions_float[index] < 0:
            positions_float[index] = positions_float[index] - (positions_float[index] / (1/scale_factor))
            new_var = Var_calculator(returns, positions_float, alpha)
            positions_float[index] = positions_float[index] + (positions_float[index] / (1/scale_factor))
            NpVaR_changes.append(new_var-prev_var)
        if positions_float[index] > 0:
            positions_float[index] = positions_float[index] + (positions_float[index] / (1/scale_factor))
            new_var = Var_calculator(returns, positions_float, alpha)
            positions_float[index] = positions_float[index] - (positions_float[index] / (1/scale_factor))
            NpVaR_changes.append(new_var-prev_var)

    return NpVaR_changes
    
'''#########################################################################################################################################################################'''





    



                   

    


































            

    












    






