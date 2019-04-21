for /l %%x in (0, 1, 9) do (
   python -O test_beasley_rfcs.py
   echo %%x
   copy test_beasley_rfcs_results.txt  test_beasley_rfcs_results_r%%x.txt
)