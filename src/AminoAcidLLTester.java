import static org.junit.jupiter.api.Assertions.*;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;

public class AminoAcidLLTester {
    String Codon = "AUUGCCUCGCUA";
    String Codon2 = "CGCUAAUUGUUC";
    String Codon3 = "CGCA";
    String Codon4 = "CGCUAAUUGUUCU";
    char[] aminoAcid1 = {'R'};
    char[] aminoAcid2 = {'E', 'G', 'V', 'A','U'};
    @Test
    public void test1() {//checks to see if the method sorts the string
        char[] aminoList = AminoAcidLL.createFromRNASequence(Codon).aminoAcidList();//given amino acid list
        assertEquals(aminoAcid2, aminoList);
    }

    @Test
    public void test2() {//checks the order
        AminoAcidLL sortList = AminoAcidLL.createFromRNASequence(Codon);
        assertTrue(AminoAcidLL.sort(sortList).isSorted());
    }

    @Test
    public void test3() {//tests to see if method prints an outprint
        System.out.print(AminoAcidLL.createFromRNASequence(Codon));
    }

    @Test
    public void test8() {//looking for output
        System.out.print(AminoAcidLL.createFromRNASequence(Codon2));
    }
    @Test
    public void test5() {//Test if codonCompare can compare the two lists
        AminoAcidLL sortList1 = AminoAcidLL.createFromRNASequence(Codon);
        AminoAcidLL sortList2 = AminoAcidLL.createFromRNASequence(Codon2);
        assertEquals(1, sortList1.codonCompare(sortList2));
    }
    @Test
    public void test6() {//checks to see if sorted
        AminoAcidLL sortList = AminoAcidLL.createFromRNASequence(Codon2);
        assertTrue(AminoAcidLL.sort(sortList).isSorted());
    }

    @Test
    public void test7() {//checks order between 4 codons
        char[] aminoList = AminoAcidLL.createFromRNASequence(Codon3).aminoAcidList();//given amino acid list
        assertEquals(aminoAcid2, aminoList);//compare expected to actual output
    }

    @Test
    public void test4() {//checks to see if contents are sorted
        AminoAcidLL sortList = AminoAcidLL.createFromRNASequence(Codon);
        System.out.print(AminoAcidLL.sort(sortList));
    }

    @Test
    public void test9() {//checks between 5 codons
        AminoAcidLL sortList = AminoAcidLL.createFromRNASequence(Codon4);
        System.out.print(AminoAcidLL.sort(sortList));
    }

    @Test
    public void test10() {//Test if codonCompare can successfully compare two lists
        AminoAcidLL sortList1 = AminoAcidLL.createFromRNASequence(Codon);
        AminoAcidLL sortList2 = AminoAcidLL.createFromRNASequence(Codon3);
        assertEquals(0, sortList1.codonCompare(sortList2));
    }
}
}
