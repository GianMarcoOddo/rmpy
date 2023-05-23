
## check_ticker_single ################################################################################################################################################
############################################################################################################################################################################

def check_ticker_single(ticker, types, name, single_value=False, values=None):
    '''
    Check if a ticker value is valid.

    Args:
    ticker: A ticker value to check.
    types: A tuple of types that the ticker value can be. 
    name: A string representing the name of the ticker value. 
    single_value: A boolean indicating whether the function should only accept a single value. 
    values: A list of valid ticker values to check against. 
    Returns:
    The validated ticker value as a string.
    Raises:
    TypeError: If the ticker value is not one of the types specified in `types`.
    ValueError: If `single_value` is `True` and more than one ticker value is provided, 
    or if `values` is provided and the ticker value is not in the list of valid values.
    '''
    import pandas as pd ; import numpy as np

    if ticker == "":
        raise TypeError("Please insert a ticker")
    
    # Check ticker type
    if not isinstance(ticker, types):
        raise TypeError(f"The {name} should be one of the following types: {types}")

    # Convert to list
    if isinstance(ticker, (list, pd.Series, np.ndarray)):
        if len(ticker) == 0:
            raise ValueError("No ticker provided")
        ticker = list(ticker)
        # Check for single value
        if single_value and len(ticker) > 1:
            raise ValueError(f"More than one {name} has been provided. This function works only with one {name} at the time")
        ticker = ticker[0]

    # Convert to upper case
    if not isinstance(ticker, (str)):
        raise TypeError("The ticker must be a string")
    else: 
        ticker = ticker.upper()

    # Check in values
    if values and ticker not in values:
        raise ValueError(f"Please, be sure to use a correct {name}! {values}")
    
    return ticker

## check_position_single ################################################################################################################################################
############################################################################################################################################################################

def check_position_single(position):
    """
    Check if a position value is valid.

    Args:
    position: A position value to check.
    types: A tuple of types that the position value can be. 
    name: A string representing the name of the position value. 
    single_value: A boolean indicating whether the function should only accept a single value. 
    values: A list of valid position values to check against. 
    Returns:
    The validated position value as a string.
    Raises:
    TypeError: If the position value is not one of the types specified in `types`.
    ValueError: If `single_value` is `True` and more than one position value is provided, 
    or if `values` is provided and the position value is not in the list of valid values.
    """
    import pandas as pd ; import numpy as np

    # Check input type
    if not isinstance(position, (int, np.int32, np.int64, float, np.float32, np.float64, list, np.ndarray, pd.Series)):
        raise TypeError(f"The position should be a number!")

    # Convert to list
    if isinstance(position, (list, pd.Series, np.ndarray)):
        position = list(position)
        if len(position) == 0: 
            raise ValueError("No position provided")
        # Check for single value
        if len(position) > 1:
            raise ValueError(f"More than one position has been provided. This function works only with one position at the time")
    
    if isinstance(position, (list)):
        position = position[0]

    if not isinstance(position, (int, np.int32, np.int64, float, np.float32, np.float64)):
        raise TypeError(f"The position should be a number!")

    return position

## check_and_convert_dates ################################################################################################################################################
############################################################################################################################################################################

def check_and_convert_dates(start_date, end_date):

    """
    Check if the input dates are valid and convert them to a datetime object.

    Args:
    start_date: A string representing the start date in one of the following formats: 
    "%Y-%m-%d", "%Y.%m.%d", "%Y/%m/%d", or "%Y_%m_%d", or the string "today" to represent today's date.
    end_date: A string representing the end date in one of the following formats: 
    "%Y-%m-%d", "%Y.%m.%d", "%Y/%m/%d", or "%Y_%m_%d", or the string "today" to represent today's date.
    Returns:
    A tuple containing the start and end dates as datetime objects.
    Raises:
    ValueError: If the start_date or end_date is not in one of the valid formats, or if the 
    end_date is before the start_date, or if the end_date is after today's date, or if the start_date is equal to the end_date.
    """
    from datetime import datetime; import dateutil.parser ; from datetime import date 

    if start_date == "today":
        raise ValueError(f"Please insert in the function a valid start_date format: [%Y-%m-%d, %Y.%m.%d, %Y/%m/%d ,%Y_%m_%d]")
    
    # Function to convert a single date
    def convert_to_date(input, name):
        if input == "today":
            return datetime.today()
        try:
            return dateutil.parser.parse(input)
        except ValueError:
            raise ValueError(f"Please insert in the function a valid {name} format: [%Y-%m-%d, %Y.%m.%d, %Y/%m/%d ,%Y_%m_%d]")

    start_date = convert_to_date(start_date, "start_date")
    end_date = convert_to_date(end_date, "end_date")

    # Get today's date
    today = datetime.today()

    if end_date < start_date:
        raise ValueError("The end_date cannot be before the start_date")
    if end_date > today:
        raise ValueError(f"The end_date cannot be after today's date {str(date.today())}")
    if start_date == end_date:
            raise ValueError("The start_date cannot be equal to the end_date")
    
    return start_date, end_date

## check_freq ################################################################################################################################################
############################################################################################################################################################################

def check_freq(freq, allowed_values=["daily", "weekly", "monthly"]):

    '''
    Validate the frequency input for a financial data analysis function.
        
    This function checks if the input frequency is of an allowed type and if it's a single value from the allowed values list.
    It raises an error if the input doesn't meet these requirements.
    Args:
    freq : str, list, np.ndarray, or pd.Series
    The input frequency value(s) to be checked.
    allowed_values : list, optional, default=["daily", "weekly", "monthly"]
    The list of allowed frequency values.
    Raises:
    TypeError:
        If the input frequency is not one of the allowed types (str, list, np.ndarray, or pd.Series).
     ValueError:
        If more than one frequency value is provided or if the input frequency value is not in the allowed values list.
        
    '''
    import pandas as pd
    import numpy as np

    def check_input_type(input, types, name):
        if not isinstance(input, types):
            raise TypeError(f"The {name} should be one of the following types: {types}")

    def check_single_value(input, name):
        if isinstance(input, (list, np.ndarray, pd.Series)) and len(input) > 1:
            raise ValueError(f"More than one {name} has been provided. This function works only with one {name} at the time")

    def convert_to_list(input):
        if isinstance(input, (list, pd.Series, np.ndarray)):
            input =  list(input)
        if isinstance(input, list):
            input =  input[0]
        return input

    def check_in_values(input, values, name):
        if not isinstance(input, str):
            raise TypeError("The freq parameter should be a string!")
        if input not in values:
            raise ValueError(f"Please, be sure to use a correct {name}! {values}")


    check_input_type(freq, (str, list, np.ndarray, pd.Series), "freq")
    check_single_value(freq, "freq")
    freq = convert_to_list(freq)
    check_in_values(freq, allowed_values, "freq")

    return freq

## check_alpha ################################################################################################################################################
############################################################################################################################################################################

def check_alpha(alpha):
    """
    Check if an alpha value or a list of alpha values is valid.
    Args:
    alpha: A float, list, numpy array, or pandas series representing one or more alpha values.
    Returns:
    alpha: A float representing a valid alpha value.
    Raises:
    TypeError: If the alpha value is not a float, list, numpy array, or pandas series.
    ValueError: If more than one alpha value is provided, or if the alpha value is not within the range of 0 to 1.
    """

    import pandas as pd
    import numpy as np

    # Check input type
    if not isinstance(alpha, (float, np.float32, np.float64, list, np.ndarray, pd.Series)):
        raise TypeError(f"The alpha should be one of the following types: float, list, np.ndarray, pd.Series")

    # Convert to list if input is pd.Series or np.ndarray
    if isinstance(alpha, (pd.Series, np.ndarray)):
        if len(alpha) == 0: 
            raise ValueError("No alpha provided")
        if len(alpha) > 1:
            raise ValueError("More than one alpha provided")
        alpha =  list(alpha)

    # Check if input is list and convert to single float
    if isinstance(alpha, list):
        if len(alpha) == 0: 
            raise ValueError("No alpha provided")
        if len(alpha) > 1:
            raise ValueError("More than one alpha provided")
        alpha = alpha[0]

    # Check alpha range
    if not isinstance(alpha, (float, np.float32, np.float64)):
        raise TypeError("The alpha value should be a number (float)")
    if alpha <= 0 or alpha >= 1:
        raise ValueError("The alpha value should be between 0 and 1")

    return alpha


## check_kind ################################################################################################################################################
############################################################################################################################################################################


def check_kind(kind, allowed_values=["abs", "rel"]):

    """
    Check if a kind value or a list of kind values is valid.
    Args:
    kind: A string, list, numpy array, or pandas series representing one or more kind values.
    allowed_values: A list of valid kind values.
    Returns:
    None
    Raises:
    TypeError: If the kind value is not a string, list, numpy array, or pandas series.
    ValueError: If more than one kind value is provided, or if the kind value is not one of the allowed values.
    """

    import pandas as pd
    import numpy as np

    if not isinstance(kind, (str, list, np.ndarray, pd.Series)):
        raise TypeError("The kind should be one of the following types: (str, list, np.ndarray, pd.Series)")

    if isinstance(kind, (pd.Series, np.ndarray)):
        kind = list(kind)

    if isinstance(kind, list):
        if len(kind) > 1:
            raise ValueError("More than one kind has been provided. This function works only with one kind at a time")
        if len(kind) == 1:
            kind = kind[0]

    if kind not in allowed_values:
        raise ValueError(f"Please, be sure to use a correct kind! {allowed_values}")
    
    return kind


## check_display ################################################################################################################################################
############################################################################################################################################################################


def check_display(display, allowed_values=[True, False]):

    """
    Check if a display value or a list of display values is valid.
    Args:
    display: A boolean, string, list, numpy array, or pandas series representing one or more display values.
    allowed_values: A list of valid display values.
    Returns:
    None
    Raises:
    TypeError: If the display value is not a boolean, string, list, numpy array, or pandas series.
    ValueError: If more than one display value is provided, or if the display value is not one of the allowed values.
    """

    import pandas as pd
    import numpy as np

    def check_input_type(input, types, name):
        if not isinstance(input, types):
            raise TypeError(f"The {name} should be one of the following types: {types}")

    def check_single_value(input, name):
        if isinstance(input, (list, np.ndarray, pd.Series)) and len(input) > 1:
            raise ValueError(f"More than one {name} has been provided. This function works only with one {name} at a time")

    def convert_to_list(input):
        if isinstance(input, (pd.Series, np.ndarray)):
            return list(input)
        if isinstance(input, list):
            return input
        return [input]

    def check_in_values(input, values, name):
        if input[0] not in values:
            raise ValueError(f"Please, be sure to use a correct {name}! {values}")

    check_input_type(display, (bool, str, list, np.ndarray, pd.Series), "display")
    display = convert_to_list(display)
    check_single_value(display, "display")
    check_in_values(display, allowed_values, "display")


## check_tickers_port ################################################################################################################################################
############################################################################################################################################################################


def check_tickers_port(tickers):
    """
    Check if a list of tickers is valid and return a list of tickers in upper case.
    Args:
    tickers: A numpy array, pandas series, or list of ticker symbols.
    Returns:
    A list of tickers in upper case.
    Raises:
    TypeError: If the tickers are provided in a pandas dataframe.
    Exception: If the tickers are an empty list, array, or series.
    Exception: If the tickers array or list is multidimensional.
    Exception: If only one ticker or a non-string type is provided.
    ValueError: If there are duplicate tickers in the list or array.
    """

    import pandas as pd ; import numpy as np

    if not isinstance(tickers, (list, np.ndarray, pd.Series)):
        raise TypeError("Tickers should be provided in the following object: [ np.ndarray, pd.core.series.Series, list ]")

    # Convert to upper case
    tickers = [x.upper() for x in tickers]

    if isinstance(tickers, pd.Series):
        tickers = list(tickers)
        if len(tickers) == 0:
            raise Exception("The Series with the tickers is empty")

    if isinstance(tickers, np.ndarray):
        if len(tickers) == 0:
            raise Exception("The array with the tickers is empty")

        if tickers.ndim > 1:
            raise Exception("The tickers array should be a one-dimensional array")
        
        tickers = list(tickers)

    if isinstance(tickers, list):

        for ticker in tickers:
            if ticker == "":
                raise TypeError("Either one ticker or more than one ticker not provided")
            
        if len(tickers) == 0:
            raise Exception("The list with the tickers is empty")
        if np.ndim(tickers) > 1:
            raise Exception("The tickers list should be a one-dimensional list")
        if len(tickers) == 1:
            raise Exception("Only one ticker has been provided, this function works with multiple tickers")
            # Check for duplicate tickers
        if len(tickers) != len(set(tickers)):
            raise ValueError("Duplicate tickers found. Please ensure all tickers are unique.")
    
    if not all(isinstance(item, str) for item in tickers):
        raise TypeError("The list, array or Series of tickers should contain only strings")

    
    return tickers


## check_positions_port ################################################################################################################################################
############################################################################################################################################################################

def check_positions_port(positions):

    """
    Check if a list of positions is valid and return a numpy array of positions.
    Args:
    positions: A numpy array, pandas series, list, or single number representing the positions.
    Returns:
    A numpy array of positions.
    Raises:
    TypeError: If the positions are provided in a pandas dataframe.
    Exception: If the positions are an empty list, array, or series.
    Exception: If the positions array or list is multidimensional.
    TypeError: If the array or list of positions contains non-numeric types.
    """

    import pandas as pd ; import numpy as np

    if not isinstance(positions, (np.ndarray, pd.core.series.Series, list)):
        raise TypeError("Positions should be provided in the following object: [ np.ndarray, pd.core.series.Series, list ]")
    
    if isinstance(positions, pd.Series):
        positions = positions.to_numpy()

        if len(positions) == 0: 
            raise Exception("The Series with the positions is empty")
        
    if isinstance(positions, np.ndarray):
        if len(positions) == 0: 
            raise Exception("The array with the positions is empty")
        
        if positions.ndim > 1:
            raise Exception("The positions array should be a one-dimensional array")
        
        if not np.all(np.isreal(positions)):
            raise TypeError("The array of positions should contain only numbers")
        
    if isinstance(positions, list):
        if len(positions) == 0: 
            raise Exception("The list with the positions is empty")
        
        if np.ndim(positions) > 1:
            raise Exception("The positions list should be a one-dimensional list")
        
        if not all(isinstance(item, (int, float)) for item in positions):
            raise TypeError("The list of positions should contain only numbers")
        
        positions = np.array(positions)

    return positions

## validate_returns_single ################################################################################################################################################
############################################################################################################################################################################

def validate_returns_single(returns):
    """
    Check if a list, pandas series, numpy array or dataframe of returns is valid and return a numpy array of returns.
    Args:
    returns: A list, pandas series, numpy array, or dataframe of returns for a single asset.
    Returns:
    A numpy array of returns.
    Raises:
    TypeError: If the returns are not provided in a list, pandas series, numpy array or dataframe.
    Exception: If the list, series, array or dataframe of returns is empty.
    Exception: If the returns list, series, array or dataframe has more than one dimension or more than one column.
    TypeError: If the returns contain non-numeric types.
    Exception: If the returns contain NaN values.
    """

    import pandas as pd ; import numpy as np

    if not isinstance(returns,(list,pd.core.series.Series,np.ndarray, pd.DataFrame)):
        raise TypeError("Returns must be provided in the form of: [list, pd.Series, np.array, pd.DataFrame]")

    if isinstance(returns,list):
        if len(returns) == 0:
            raise Exception("The list of returns is empty")
        dim = np.ndim(returns)
        if dim > 1:
            raise Exception("The returns list should be one-dimensional or should contain returns for only one asset")
        if not all(isinstance(item, (float, np.float16, np.float32, np.float64)) for item in returns):
            raise TypeError("The list of returns should contain only numbers - (Percentages)")
        returns = np.asarray(returns)

    if isinstance(returns,pd.core.series.Series):
        returns = returns.values
        if len(returns) == 0:
            raise Exception("The Series of returns is empty")
        if not all(isinstance(item, (float, np.float16, np.float32, np.float64)) for item in returns):
            raise TypeError("The Series of returns should contain only numbers - (Percentages)")
        returns = np.asarray(returns)

    if isinstance(returns,np.ndarray):
        if len(returns) == 0:
            raise Exception("The array of returns is empty")
        dim = np.ndim(returns)
        if dim > 1:
            raise Exception("The returns array should be one-dimensional or should contain returns for only one asset")
        if not all(isinstance(item, (float, np.float16, np.float32, np.float64)) for item in returns):
            raise TypeError("The array of returns should contain only numbers")

    if isinstance(returns,pd.DataFrame):
        if returns.empty:
            raise Exception("The DataFrame with the returns is empty")
        if returns.shape[1] > 1:
            raise Exception("A DataFrame with more than one column has been provided, this function works only with one asset at the time")
        for col in returns.columns:
            if not all(isinstance(item, (float, np.float16, np.float32, np.float64)) for item in returns[col]):
                raise TypeError(f"The DataFrame column {col} should contain only numbers - Percentages")
        if returns.shape[1] == 1:
            returns = returns.values
    
    # After converting to numpy array, check for NaN values
    if np.isnan(returns).any():
        raise Exception("The returns contain NaN values")

    return returns

## validate_returns_port ################################################################################################################################################
############################################################################################################################################################################

def validate_returns_port(returns):

    """
    Validates the format of the returns for a portfolio of assets.
    Args:
    returns : np.ndarray or pd.DataFrame with the returns of a portfolio of assets. 
    Must be provided in the form of an np.array or pd.DataFrame.
    Returns:
    np.ndarray with the returns of a portfolio of assets as a np.ndarray.
    Raises:
    Exception
        If the returns are not provided as an np.array or pd.DataFrame.
        If the np.array with the returns contains less than one column. 
        If the np.array with the returns contains NaN values.
        If the DataFrame with the returns is empty.
        If the DataFrame with the returns contains less than two columns.
        If the DataFrame with the returns contains NaN values.
        If any column of the DataFrame with the returns contains values that are not numbers.
    """

    import pandas as pd ; import numpy as np

    if not isinstance(returns,(np.ndarray, pd.DataFrame)):
        raise Exception("Returns must be provided in the form of: [np.array, pd.DataFrame]. This function works with multiple assets")

    # Array
    if isinstance(returns,np.ndarray):
        if returns.shape[1] < 1:
            raise Exception("The array with the returns contains less than or one column. This function works with more than one asset")

        if np.isnan(returns).any():
            raise Exception("The array with the returns contains NaN values")

        for col in range(returns.shape[1]):
            if not all(isinstance(item, (int, float)) for item in returns[:,col]):
                raise Exception("The Array of returns should contain only numbers")

    # DataFrame
    if isinstance(returns,pd.DataFrame):
        if returns.empty:
            raise Exception("The DataFrame with the returns is empty")

        if returns.shape[1] == 1:
            raise Exception("The DataFrame with the returns should have more than one column")

        if returns.isna().values.any():
            raise Exception("The DataFrame with the returns contains NaN values")

        for col in returns.columns:
            if not all(isinstance(item, (int, float)) for item in returns[col]):
                raise Exception(f"The DataFrame column {col} should contain only numbers")

        returns = returns.values

    return returns


## check_scale_factor ################################################################################################################################################
############################################################################################################################################################################

def check_scale_factor(scale_factor):

    """
    Check if an scale_factor value or a list of scale_factor values is valid.
    Args:
    scale_factor: A float, list, numpy array, or pandas series representing one or more scale_factor values.
    Returns:
    scale_factor
    Raises:
    TypeError: If the scale_factor value is not a float, list, numpy array, or pandas series.
    ValueError: If more than one scale_factor value is provided, or if the scale_factor value is not within the range of 0 to 1.
    """

    import pandas as pd
    import numpy as np

    def check_input_type(input, types, name):
        if not isinstance(input, types):
            raise TypeError(f"The {name} should be one of the following types: {types}")

    def check_single_value(input, name):
        if isinstance(input, (list, np.ndarray, pd.Series)) and len(input) > 1:
            raise ValueError(f"More than one {name} has been provided. This function works only with one {name} at a time")

    def convert_to_list(input):
        if isinstance(input, (pd.Series, np.ndarray)):
            if len(input) == 0: 
                raise ValueError("No scale_factor provided")
            if len(input) > 1:
                raise Exception("More than one scale_factor provided")
            input=  list(input)
        if isinstance(input, list):
            if len(input) == 0: 
                raise ValueError("No scale_factor provided")
            if len(input) > 1:
                raise Exception("More than one scale_factor provided")
            input =  input[0]
        return input

    def check_scale_factor_range(scale_factor):
        if not isinstance(scale_factor, (int, float)):
            raise TypeError("The scale_factor value should be a number (float)")
        if scale_factor <= 0 or scale_factor >= 1:
            raise ValueError("The scale_factor value should be between 0 and 1")

    check_input_type(scale_factor, (float, list, np.ndarray, pd.Series), "scale_factor")
    scale_factor = convert_to_list(scale_factor)
    check_single_value(scale_factor, "scale_factor")
    check_scale_factor_range(scale_factor)

    return scale_factor


## check_confidence_level ################################################################################################################################################
############################################################################################################################################################################

def check_confidence_level(confidence_level):
    """
    This function checks if the provided confidence level parameter is valid for statistical analysis.
    Args:
    confidence_level: (np.ndarray, pd.core.series.Series, list, number)
        The confidence level parameter to be checked for validity.
    Returns:
    confidence_level: (float)
        The validated confidence level to be used for statistical analysis.
    Raises:
    TypeError: 
        - If the confidence_level parameter is a pandas DataFrame object.
        - If the confidence_level parameter is not a float number or an array with all float values.
    Exception:
        - If the Series or array/list with the confidence_level parameter is empty.
        - If the confidence_level array/list is more than one-dimensional.
        - If there is more than one confidence_level parameter provided.
        - If the confidence_level parameter is not a float number between 0 and 1.
    """

    import pandas as pd ; import numpy as np

    # Raise TypeError if confidence_level is a DataFrame
    if isinstance(confidence_level, pd.core.frame.DataFrame):
        raise TypeError("The confidence_level parameter should be provided in the following object: [np.ndarray, pd.core.series.Series, list, number]")

    # Handle if confidence_level is a Series
    if isinstance(confidence_level, pd.core.series.Series):
        if len(confidence_level) == 0:
            raise Exception("The Series with the confidence_level parameter is empty")
        confidence_level = list(confidence_level)

    # Handle if confidence_level is a list or ndarray

    if isinstance(confidence_level, (np.ndarray, list)):
        if len(confidence_level) == 0:
            raise Exception("The array/list with the confidence_level parameter is empty")
        if len(confidence_level) > 1:
            raise Exception("The array/list with the confidence_level parameter should contain only one value")
        dim = np.ndim(confidence_level)
        if dim > 1:
            raise Exception("The confidence_level array/list should be one-dimensional")
        confidence_level = confidence_level[0]

    # Handle if confidence_level is a number
    if not isinstance(confidence_level, (float, np.int32, np.int64, np.float16, np.float32, np.float64)):
        raise TypeError("The confidence_level parameter must be a float number between 0 and 1")

    # Ensure the value of confidence_level is between 0 and 1
    if confidence_level > 1 or confidence_level < 0:
        raise Exception("Please, insert a correct value for confidence_level parameter! Between 0 and 1 (i.e for 95% insert 0.95)")

    return confidence_level

## check_interval ################################################################################################################################################
############################################################################################################################################################################

def check_interval(interval):
    """
    This function checks if the provided interval parameter is valid for statistical analysis.
    Args:
    interval: (np.ndarray, pd.core.series.Series, list, number)
        The interval parameter to be checked for validity.
    Returns:
    interval: (int)
        The validated interval to be used for statistical analysis.
    Raises:
    TypeError: 
        - If the interval parameter is a pandas DataFrame object.
        - If the interval parameter is not an integer or an array with all integer values.
    Exception:
        - If the Series or array/list with the interval parameter is empty.
        - If the interval array/list is more than one-dimensional.
        - If there is more than one interval parameter provided.
        - If the interval parameter is not an integer greater than 1.
    """
    import pandas as pd ; import numpy as np

    # Raise TypeError if interval is a DataFrame
    if isinstance(interval, pd.core.frame.DataFrame):
        raise TypeError("The interval parameter should be provided in the following object: [np.ndarray, pd.core.series.Series, list, number]")

    # Handle if interval is a Series
    if isinstance(interval, pd.core.series.Series):
        if len(interval) == 0:
            raise Exception("The Series with the interval parameter is empty")
        interval = list(interval)

    # Handle if interval is a list or ndarray
    if isinstance(interval, (np.ndarray, list)):
        if len(interval) == 0:
            raise Exception("The array/list with the interval parameter is empty")
        dim = np.ndim(interval)
        if dim > 1:
            raise Exception("The interval array/list should be one-dimensional")
        if len(interval) > 1:
            raise Exception("More than one interva provided")
        interval = interval[0]
        

    # Handle if interval is a number
    if not isinstance(interval, (int, np.int32, np.int64 )):
        raise TypeError("The interval parameter must be an integer ")

    # Ensure the value of interval higher than one
    if interval < 1:
        raise Exception("Please, insert a correct value for interval parameter! > 1)")

    return interval

""" ########################################################################################################################### """

