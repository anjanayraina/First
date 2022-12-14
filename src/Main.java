import java.io.*;
import java.util.*;

public class Main{
    static int MAX_SEQUENCE_LENGTH = 10;
    static int  PSEUDOCOUNT = 2;
   static  class Wrapper{
        double[][] profile;
        char[] amino_acids;
        Wrapper(double[][] profile , char[] amino_acids){
            this.profile = profile;
            this.amino_acids = amino_acids;
        }
    }

    public static char[] find_amino_acids(String sequences[]) {
        ArrayList<Character> amino_acids = new ArrayList<>();
        for (String seq : sequences) {
            for(char c : seq.toCharArray()){
                if(!amino_acids.contains(c)){
                    amino_acids.add( c);
                }
            }
        }
        int  n =amino_acids.size();
        char[] res=  new char[n];
        for(int i=0;i<n;i++){
            res[i] = amino_acids.get(i);
        }
//        for(char c : res) System.out.print(c + " ");
        return res;
    }

    public static int find_num_occurrences(char amino, int index, String sequences[]) {
        int counter = 0;
        for (String seq : sequences) {
            if (amino == seq.charAt(index)) {
                counter++;
            }
        }
        return counter;
    }
    public static Wrapper create_PSSM_matrix(String sequences[]){
        int num_of_seqs = sequences.length;
        int seq_length = sequences[0].length();
        char[] amino_acids = find_amino_acids(sequences);
        double profile[][] = new double[amino_acids.length][seq_length];
        for(int  amino_i=0;amino_i<amino_acids.length;amino_i++){
            for(int i=0;i<seq_length;i++){
                profile[amino_i][i] = find_num_occurrences(amino_acids[amino_i] , i , sequences);
            }
        }

        for(int i=0;i< profile.length;i++){
            for(int j =0;j<profile[0].length;j++){
                profile[i][j] +=PSEUDOCOUNT;
            }
        }


        for(int i=0;i< profile.length;i++){
            for(int j =0;j<profile[0].length;j++){
                profile[i][j] =profile[i][j]/(num_of_seqs + amino_acids.length * PSEUDOCOUNT);
            }
        }

        for(int row_i=0;row_i<profile.length;row_i++){
            double sum=0;
            for(int i=0;i<profile[row_i].length;i++){
                sum+=profile[row_i][i];
            }
            for(int i=0;i<profile[row_i].length;i++){
                profile[row_i][i] = profile[row_i][i] /(sum / seq_length);
            }

        }


        for(int i=0;i< profile.length;i++){
            for(int j =0;j<profile[0].length;j++){
                profile[i][j] =Math.log(profile[i][j]);
            }
        }


        return new Wrapper(profile, amino_acids);
    }

    public static int getOccurence(char []aminoAcids , char c){
       int i=0;
       for( ; i <aminoAcids.length;i++){
           if(aminoAcids[i]  == c)return i;
       }
       return -1;

    }
    public static int calculate_score(String sequence, double [][] profile, char[] aminoAcids) {
        int score = 0;

        for(int char_i =0;char_i < sequence.length();char_i++){
            int amino_index =getOccurence(aminoAcids  , sequence.charAt(char_i));
            score += amino_index;
        }

        return score;
    }

    public static List<String> add_one_gap(String [] sequences) {
        List<String> newSeqs = new ArrayList<>();
        for (String seq : sequences) {
            for (int j = 0; j < seq.length(); j++) {
                String newSeq = seq.substring(0, j) + "-" + seq.substring(j);
                if (!newSeqs.contains(newSeq)) {
                    newSeqs.add(newSeq);
                }
            }

            String newSeq = seq + "-";
            if (!newSeqs.contains(newSeq)) {
                newSeqs.add(newSeq);
            }
        }

        return newSeqs;
    }

    public static List<String> generate_gaps(String sequence, int numGaps) {
        List<String> finalSeqs = new ArrayList<>();
        List<String> sequences = new ArrayList<>();
        sequences.add(sequence);
        // sequences = [sequence]
        while (true) {
            String temp[] = new String[]{sequences.get(0)};
            List<String> generatedSeqs = add_one_gap(temp);
            if (generatedSeqs.get(0).length() == sequence.length() + numGaps) {
                finalSeqs.addAll(generatedSeqs);
                break;
            }
            sequences = generatedSeqs;
        }

        return finalSeqs;
    }

    public static List<String> create_seqs_with_gap(String sequence, int max_length){
        int max_num_of_gaps = max_length - sequence.length();
        List<String> generated_sequences = generate_gaps(sequence, max_num_of_gaps);
        return generated_sequences;
    }

    public static String find_best_subsequence(double[][] profile, String sequence, char[] aminoAcids){
        float maxScore = Float.MIN_VALUE;
        String bestSubseq = "";
        int seqLength = profile[1].length;
        for(int sl = seqLength ; sl > 0 ; sl--){
            for(int i = 0 ; i < sequence.length() - sl + 1 ; i++){
                String tmpSeq = sequence.substring(i, i + sl  );
                int score = calculate_score(tmpSeq, profile, aminoAcids);
                if(score > maxScore){
                    bestSubseq = tmpSeq;
                    maxScore = score;
//                    System.out.println("if entered");
                }

                // here I create different versions with gap
                if(tmpSeq.length() < seqLength){
                    List<String> tmpsWithGap = new ArrayList<String>();
                    tmpsWithGap = create_seqs_with_gap(tmpSeq, seqLength);

                    for(String seq : tmpsWithGap){
                        int scoreWithGap = calculate_score(seq, profile, aminoAcids);
                        if (scoreWithGap > maxScore){
                            bestSubseq = seq;
                            maxScore = scoreWithGap;
                        }
                    }
                }
            }
        }
        System.out.println(bestSubseq);

        return bestSubseq;
    }
    public static void main(String[] args) throws Exception{
        System.out.println("Add target sequence at last in the file");
        ArrayList<String> sequences = new ArrayList<String>();
        File file = new File("input.txt");
        BufferedReader br = new BufferedReader(new FileReader(file));

        String sequence;
        while((sequence = br.readLine()) != null){
            sequences.add(sequence);
            // System.out.println(sequence);
        }
//        Scanner scn = new Scanner(System.in);
//        int n = scn.nextInt();
//        for(int  i=0;i<n;i++){
//            sequences.add(scn.next());
//        }

        String targetSequence = sequences.get(sequences.size()-1);
        sequences.remove(sequences.size()-1);
        // all the sequences are in sequences in arraylist
        // now converting arraylsit to array
        String[] arr = new String[sequences.size()];
        arr = sequences.toArray(arr);
        System.out.println("Sequences are:\n**************");
        for(String s : arr){
            System.out.println(s);
        }
        System.out.println("**************");
        // now all the sequences are in arr array
        System.out.println("Target Sequence: " + targetSequence);

        Wrapper w = create_PSSM_matrix(arr);

        String best_subsequence =find_best_subsequence(w.profile , targetSequence, w.amino_acids);
        System.out.println("The best Sequence is : " + best_subsequence);
        System.out.println();
    }
}