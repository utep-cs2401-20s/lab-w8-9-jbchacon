class AminoAcidLL{
  public static AminoAcidLL tempHead;
  public static AminoAcidLL head;
  private char aminoAcid;
  private String[] codons;
  private int[] counts;

  private AminoAcidLL next;

  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    incrementCodons(inCodon);
    next = null;

  }
  private void incrementCodons(String inCodon){
    for(int i = 0; i<codons.length; i++){
      if(codons[i].equals(inCodon)){
        counts[i]++;// increments counter
      }
    }

  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops,
   * if not passes the task to the next node.
   * If there is no next node, add a new node to the list that would contain the codon.
   */
  private void addCodon(String inCodon) {

    if (aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
      incrementCodons(inCodon);//increment counter
    } else {
      if (next != null) {
        next.addCodon(inCodon);//new node
      } else {
        AminoAcidLL newAmino = new AminoAcidLL(inCodon);
        next = newAmino;
      }
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int total = 0;
    for(int i = 0; i<counts.length; i++){
      total+= counts[i];//adds all elements
    }
    return total;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts.
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    if(inList.isSorted()) {
      if (next == null && inList.next == null)
        return totalCount() - totalDiff(inList);

      if (inList.next == null)
        return this.totalCount() + next.aminoAcidCompare(inList);

      if (next == null)
        return inList.totalCount() + this.totalDiff(inList.next);
    }
    return this.totalDiff(inList) + next.aminoAcidCompare(inList.next);

  }

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    if(inList.isSorted()) {
      if (next == null && inList.next == null)
        return codonDiff(inList);

      if (inList.next == null)
        return this.totalCount() + next.codonCompare(inList);

      if (next == null)
        return inList.totalCount() + this.codonDiff(inList.next);
    }
    return this.totalDiff(inList) + next.codonCompare(inList.next);

  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    if(next == null)
      return new char[] {aminoAcid};
    else {
      char[] node = next.aminoAcidList();
      char[] list = new char[node.length + 1];//tail
      //for loop to fill new char array
      for (int i = 0; i < node.length; i++) {
        list[0] = aminoAcid;// head
        list[i + 1] = node[i];//every element
      }
      return list;
    }
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts() {
    if (next == null)
      return new int[]{};
    else {
      int[] node = next.aminoAcidCounts();
      int[] counts = new int[node.length + 1];

      //for loop to fill new counts array
      for (int i = 0; i < node.length; i++) {
        counts[0] = aminoAcid;// head
        counts[i + 1] = node[i];//every element
      }
      return aminoAcidCounts();
    }
  }

  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted() {
    if (tempHead == null || tempHead.next == null)
      return true;

    if (tempHead.aminoAcid > next.aminoAcid){
      return false;
  }
  tempHead = tempHead.next;
     next.isSorted();//recursion
      return true;
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence) {
    String codon = inSequence.substring(0, 3);// Every 3 characters in codon
    boolean test = true;// checks for a stop
    tempHead = head;//creates head
    if (inSequence.substring(0, 3).charAt(0) == '#') {
      head.addCodon(inSequence.substring(0, 3));
      test = false;
    } else {
      head.addCodon(inSequence.substring(0, 3));
    }
    for (int i = 3; i < inSequence.length() - 2; i += 3) {
      if (inSequence.charAt(1) == '#') {
        head.addCodon(inSequence.substring(i, i + 3));
        test = false;
      } else {
        head.addCodon(inSequence.substring(i, i + 3));
      }
    }
    return head;
  }
  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){

    if(inList.isSorted())
      return inList;
    else{
      char temp;
      for(AminoAcidLL i = inList; i.next != null; i = i.next) {
        for (AminoAcidLL j = i.next; j.next != null; j = j.next) {
          if (i.aminoAcid > j.aminoAcid) {
            temp = i.aminoAcid;
            j.next = i;
            i.aminoAcid = temp;
          }
        }
      }
    return null;// sorted
  }
 }
}
