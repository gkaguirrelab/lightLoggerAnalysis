def ndf_to_lux(cal_filepath: str,
               force_recalc: bool = False,
               matlab_engine: object | None = None
            ) -> float | int:

    # Import the MATLAB engine if not passed in and initialize it.
    need_to_initialize_engine: bool = matlab_engine is None
    if(need_to_initialize_engine):
        import matlab.engine

        matlab_engine = matlab.engine.start_matlab()
        matlab_engine.tbUseProject("lightLoggerAnalysis", nargout=0)

    try:
        # Call the helper function to calculate the lux.
        lux: float | int = matlab_engine.ndf2lux(
            cal_filepath,
            "force_recalc",
            force_recalc,
            nargout=1,
        )
    finally:
        # Close the engine if we defined it here.
        if(need_to_initialize_engine):
            matlab_engine.quit()

    return lux


def main():
    pass 

if(__name__ == "__main__"):
    main() 
