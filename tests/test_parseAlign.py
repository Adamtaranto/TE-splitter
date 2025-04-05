import os
import tempfile

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tsplit.parseAlign import filterCoordsFileTIR


def test_filterCoordsFileTIR():
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tempdir:
        # Create a temporary coords file
        coords_file_path = os.path.join(tempdir, 'TIR_element_blastn.coords')
        with open(coords_file_path, 'w') as f:
            f.write(
                '1\t4963\t1\t4963\t4963\t4963\t100.000\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '769\t1336\t644\t1210\t568\t512\t90.141\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '644\t1210\t769\t1336\t568\t512\t90.141\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '889\t1355\t637\t1103\t468\t408\t87.179\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '637\t1103\t889\t1355\t468\t408\t87.179\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '1015\t1355\t637\t977\t342\t299\t87.427\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '637\t977\t1015\t1355\t342\t299\t87.427\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '1141\t1355\t637\t851\t215\t184\t85.581\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '637\t851\t1141\t1355\t215\t184\t85.581\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '4868\t4963\t96\t1\t96\t78\t81.250\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '1\t106\t4963\t4860\t106\t85\t80.189\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '4520\t4592\t501\t429\t73\t59\t80.822\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '429\t504\t4592\t4517\t76\t61\t80.263\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '335\t374\t4684\t4645\t40\t32\t80.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '4739\t4752\t4727\t4714\t14\t14\t100.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '4714\t4727\t4752\t4739\t14\t14\t100.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '1694\t1709\t1571\t1585\t16\t15\t93.750\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '1571\t1585\t1694\t1709\t16\t15\t93.750\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '3850\t3864\t1961\t1976\t16\t15\t93.750\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '3714\t3725\t2130\t2119\t12\t12\t100.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '3516\t3530\t2709\t2695\t15\t14\t93.333\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '2695\t2709\t3530\t3516\t15\t14\t93.333\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '2119\t2130\t3725\t3714\t12\t12\t100.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '3719\t3739\t3739\t3719\t21\t18\t85.714\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '1961\t1976\t3850\t3864\t16\t15\t93.750\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '4530\t4544\t1387\t1374\t15\t14\t93.333\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '3467\t3477\t1634\t1624\t11\t11\t100.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '3969\t3985\t2019\t2003\t17\t15\t88.235\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '4255\t4265\t2265\t2255\t11\t11\t100.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '3198\t3208\t2612\t2622\t11\t11\t100.000\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '4735\t4745\t2642\t2652\t11\t11\t100.000\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '2933\t2943\t2716\t2726\t11\t11\t100.000\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '2716\t2726\t2933\t2943\t11\t11\t100.000\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '3525\t3538\t2974\t2987\t14\t13\t92.857\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '2612\t2622\t3198\t3208\t11\t11\t100.000\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '1624\t1634\t3477\t3467\t11\t11\t100.000\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '2974\t2987\t3525\t3538\t14\t13\t92.857\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
                '2006\t2019\t3982\t3969\t14\t13\t92.857\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '2255\t2270\t4265\t4248\t18\t16\t88.889\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '1374\t1387\t4544\t4530\t15\t14\t93.333\t4963\t4963\t1\t-1\tTIR_element\tTIR_element\n'
                '2642\t2652\t4735\t4745\t11\t11\t100.000\t4963\t4963\t1\t1\tTIR_element\tTIR_element\n'
            )

        # Create a mock SeqRecord
        record = SeqRecord(Seq('A' * 4963), id='TIR_element')

        # Call the function
        alignments = filterCoordsFileTIR(
            coords_file_path, record, minterm=10, flankdist=10
        )

        # Check the results
        assert len(alignments) == 1
        assert alignments[0].ref_start == 0
        assert alignments[0].ref_end == 105
        assert alignments[0].qry_start == 4962
        assert alignments[0].qry_end == 4859
        assert alignments[0].hit_length_ref == 106
        assert alignments[0].hit_length_qry == 85
        assert alignments[0].percent_identity == 80.189
        assert alignments[0].ref_name == 'TIR_element'
        assert alignments[0].qry_name == 'TIR_element'
