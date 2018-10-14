
./iDist -C ../test/2D/testHP.conf results/testHP/ -B ../test/2D/test_2_10.txt -QN ../test/2D/test_2_10_qn.txt 10

./iDist -C ../test/2D/testRefs.conf ../test/2D/test_2_10_refs.txt results/testRefs/ -B ../test/2D/test_2_10.txt -QN ../test/2D/test_2_10_qn.txt 10

./iDist -C ../test/2D/testAssigns.conf ../test/2D/test_2_10_refs.txt ../test/2D/test_2_10_assigns.txt  results/testAssigns/ -B ../test/2D/test_2_10.txt -QN ../test/2D/test_2_10_qn.txt 10

