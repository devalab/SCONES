# https://stackoverflow.com/a/66195001
def any_of(*expected_conditions):
    """ An expectation that any of multiple expected conditions is true.
    Equivalent to a logical 'OR'.
    Returns results of the first matching condition, or False if none do. """
    def any_of_condition(driver):
        for expected_condition in expected_conditions:
            from selenium.common.exceptions import WebDriverException
            try:
                result = expected_condition(driver)
                if result:
                    return result
            except WebDriverException:
                pass
        return False
    return any_of_condition