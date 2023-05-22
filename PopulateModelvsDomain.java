import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.util.Scanner;
import java.util.Arrays;
import java.util.StringTokenizer;
import javafx.util.Pair;

import java.util.regex.*;

import java.sql.*;
import gov.lbl.scop.local.LocalSQL;
import gov.lbl.scop.util.RAF;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.*;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.ce.CeParameters;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.jama.Matrix;

import java.util.concurrent.*;

public class PopulateModelvsDomain {
	
    //number of threads used in parallelization
    static int NUMTHREADS = 28;
    final public static long TIMEOUT = 600L;

    public static void splitHitTest() {
	    int[] aln0 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52};
	    int[] aln1 = {65, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 140, 141, 142, 143, 144, 145, 146, 147, 148};
	    int[][] optAlnBlock = new int[][]{aln0, aln1};
	    int[][][] optAln = new int[][][]{optAlnBlock};
	    ArrayList<ArrayList<int[]>> hits = splitHits(optAln, 15); 
	    System.out.println("hits: " + hits.get(0).get(0)[0]);
	    System.out.println("hits: " + hits.get(0).get(0)[1]);
	    System.out.println("hits: " + hits.get(1).get(0)[0]);
	    System.out.println("hits: " + hits.get(1).get(0)[1]);
	    
	    System.out.println("hits: " + hits.get(0).get(1)[0]);
	    System.out.println("hits: " + hits.get(0).get(1)[1]);
	    System.out.println("hits: " + hits.get(1).get(1)[0]);
	    System.out.println("hits: " + hits.get(1).get(1)[1]);

    }

    public static void shiftedTest() throws Exception {
	    AFPChain comp = getComparison("a.pdb", "b.pdb");
	    int[][][] optAln = comp.getOptAln();
	    int[][] optAlnBlock = optAln[0];
	    int[] aln0 = optAlnBlock[0]; // domain ?
	    int[] aln1 = optAlnBlock[1]; // model ?
	    System.out.println("aln0"+ Arrays.toString(aln0));
	    System.out.println("aln1"+ Arrays.toString(aln0));
	}


    /**
    ce
    **/
    public static AFPChain getComparison(String modelDir, String domainDir) throws Exception {

        PDBFileReader pdbreader = new PDBFileReader();
        Structure domainStructure = pdbreader.getStructure(domainDir); 
        Structure modelStructure = pdbreader.getStructure(modelDir); 

        // Fetch CA atoms for the structures to be aligned
        Atom[] caDomain = StructureTools.getAtomCAArray(domainStructure);
        Atom[] caModel = StructureTools.getAtomCAArray(modelStructure);

        // Get StructureAlignment instance
        StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
        // get default parameters
        CeParameters params = new CeParameters();
        // Perform the alignment
        return algorithm.align(caDomain, caModel, params);

    }


    /**
       Translate an index (0-based) in ATOM records of a domain seq into an index in the RAF (0-based).  -1 if not
       found.
	nodeID: scop_node.id of the domain
	ATOMDomainIdx: 0-indexed into the ATOM records of the domain sequence
	*/

    private static ArrayList<Pair<Integer, Integer>> translateATOMDomainIdx(int nodeID, int[] ATOMDomainIndices) {
	Pair<SequenceFragmentRAFIdx, SequenceFragmentRAFIdx> domainFragPair = domainSeq(nodeID, 1, 0);
	SequenceFragmentRAFIdx domainFragATOM = domainFragPair.getKey();
	SequenceFragmentRAFIdx domainFragSEQRES = domainFragPair.getValue();
	ArrayList<Pair<Integer, Integer>> rv = new ArrayList<>();
	for (Integer ATOMDomainIdx : ATOMDomainIndices) {
	    Integer RAFidx = -1;
	    Integer SEQRESidx = -1;
	    try {
	        RAFidx = domainFragATOM.RAFindices.get(ATOMDomainIdx);
	        for (int i = 0; i < domainFragSEQRES.RAFindices.size(); i++) {
		    if (domainFragSEQRES.RAFindices.get(i).equals(RAFidx)) {
			SEQRESidx = i;
			break;
		    }
	        }
	    } finally {
	        rv.add(new Pair<> (RAFidx, SEQRESidx));
	    }
	}
	return rv; 
    }

    // first axis: which sequence
    // second axis: which consecutive run
    // third axis: start and end
    private static ArrayList<ArrayList<int[]>> splitHits(int[][][] optAln, int gapTolerance) {
	    int[][] optAlnBlock = optAln[0];
	    int[] aln0 = optAlnBlock[0]; // domain ?
	    int[] aln1 = optAlnBlock[1]; // model ?
	    // System.out.println(Arrays.toString(aln0));
	    // System.out.println(Arrays.toString(aln1));
	    ArrayList<int[]> hits0 = new ArrayList<int[]>();
	    ArrayList<int[]> hits1 = new ArrayList<int[]>();
	    int start0 = aln0[0];
	    int start1 = aln1[0];
	    for (int i = 0; i < aln0.length - 1; i++) {
		    if (i == aln0.length - 2) {
			   int end0 = aln0[i] + 1; // inclusive
			   int end1 = aln1[i] + 1; // inclusive
			   hits0.add(new int[]{start0, end0});
			   hits1.add(new int[]{start1, end1});
		    } else if (!(aln1[i + 1] - aln1[i] < gapTolerance)) {
			   int end0 = aln0[i]; // inclusive
			   int end1 = aln1[i]; // inclusive
			   hits0.add(new int[]{start0, end0});
			   hits1.add(new int[]{start1, end1});
			   start0 = aln0[i + 1];
			   start1 = aln1[i + 1];
		    }
	    }
	    ArrayList<ArrayList<int[]>> hits =  new ArrayList<>();
	    hits.add(hits0);
	    hits.add(hits1);
	    return hits;
    }

    public static void populateOneComparison(String modelDir, int testModelId, String sid, int EQRThreshold, int releaseID) throws Exception {
	    	    System.out.println("model: " + testModelId + " sid: " + sid);
		    PreparedStatement stmt1 = LocalSQL.prepareStatement("select astral_domain.id as ad_id, scop_node.id as node_id from astral_domain JOIN scop_node ON astral_domain.node_id = scop_node.id JOIN link_pdb on scop_node.id = link_pdb.node_id where scop_node.sid = ? and scop_node.release_id = ? and astral_domain.style_id = 1 and astral_domain.source_id = 1;");
		    PreparedStatement stmt3 = LocalSQL.prepareStatement("insert into model_vs_domain_structure_alignment (model_id, domain_id, structure_aligner_id, z_score, model_start, model_end, domain_start_ATOM, domain_end_ATOM, domain_start_raf, domain_end_raf, domain_start_SEQRES, domain_end_SEQRES, translate_x, translate_y, translate_z, rotate_1_1, rotate_1_2, rotate_1_3, rotate_2_1, rotate_2_2, rotate_2_3, rotate_3_1, rotate_3_2, rotate_3_3, num_eq_res) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);");
			    stmt1.setString(1, sid);
			    stmt1.setInt(2, releaseID);
			    System.out.println(stmt1);
			    ResultSet domainRs = stmt1.executeQuery();
			    if (!domainRs.next()) {
				    System.out.println("No matching SCOP node or RAFv3 line for this sid");
				    return;
			    }
			    int astralDomainID = domainRs.getInt("ad_id");
			    int nodeID = domainRs.getInt("node_id");
			    if (domainRs.next()) {
				    System.out.println("There are multiple nodeIDs!");
			    }
			    String pdbStyleDir = String.join("/", "/lab/proj/astral/pdbstyle/2.06", sid.substring(2,4), sid) + ".ent";
			    AFPChain comp = getComparison(modelDir, pdbStyleDir);
			    int NrEQR = comp.getNrEQR();

			    if (NrEQR < EQRThreshold) {
				    return;
			    }

			    Matrix rotMat = comp.getBlockRotationMatrix()[0];
			    Atom shiftVec = comp.getBlockShiftVector()[0];
			    double zScore = comp.getProbability();
			    int[][][] optAln = comp.getOptAln();
			    ArrayList<ArrayList<int[]>> hits = splitHits(optAln, 50);
			    // System.out.println("hits: " + hits.get(0).get(0)[0]);
			    // System.out.println("hits: " + hits.get(0).get(0)[1]);
			    // System.out.println("hits: " + hits.get(1).get(0)[0]);
			    // System.out.println("hits: " + hits.get(1).get(0)[1]);
			    for (int h = 0; h < hits.get(0).size(); h++) {
				    int modelStart = hits.get(1).get(h)[0]; // check whether model or hit is first
				    int modelEnd = hits.get(1).get(h)[1];
				    // System.out.println("optAln: " + Arrays.deepToString(optAln));
				    // System.out.println("hit index: " + hits.get(1).get(h)[1]);
				    // System.out.println("modelEnd: " + modelEnd);
				    int[] domainBoundariesATOM = hits.get(0).get(h);
				    ArrayList<Pair<Integer,Integer>> domainBoundariesPair = translateATOMDomainIdx(nodeID, domainBoundariesATOM);
				    int domainStartRAF = domainBoundariesPair.get(0).getKey();
				    int domainEndRAF = domainBoundariesPair.get(1).getKey();
				    int domainStartSEQRES = domainBoundariesPair.get(0).getValue();
				    int domainEndSEQRES = domainBoundariesPair.get(1).getValue();
				    
				    int i = 1;
				    stmt3.setInt(i++, testModelId);
				    stmt3.setInt(i++, astralDomainID);
				    stmt3.setInt(i++, 1);
				    stmt3.setDouble(i++, zScore);
				    stmt3.setInt(i++, modelStart);
				    stmt3.setInt(i++, modelEnd);
				    stmt3.setInt(i++, domainBoundariesATOM[0]);
				    stmt3.setInt(i++, domainBoundariesATOM[1]);
				    stmt3.setInt(i++, domainStartRAF);
				    stmt3.setInt(i++, domainEndRAF);
				    stmt3.setInt(i++, domainStartSEQRES);
				    stmt3.setInt(i++, domainEndSEQRES);
				    stmt3.setDouble(i++, shiftVec.getX());
				    stmt3.setDouble(i++, shiftVec.getY());
				    stmt3.setDouble(i++, shiftVec.getZ());
				    stmt3.setDouble(i++, rotMat.get(0, 0));
				    stmt3.setDouble(i++, rotMat.get(0, 1));
				    stmt3.setDouble(i++, rotMat.get(0, 2));
				    stmt3.setDouble(i++, rotMat.get(1, 0));
				    stmt3.setDouble(i++, rotMat.get(1, 1));
				    stmt3.setDouble(i++, rotMat.get(1, 2));
				    stmt3.setDouble(i++, rotMat.get(2, 0));
				    stmt3.setDouble(i++, rotMat.get(2, 1));
				    stmt3.setDouble(i++, rotMat.get(2, 2));
				    stmt3.setDouble(i++, NrEQR);
				    stmt3.executeUpdate();
			    }
    }

	public static class SequenceFragmentRAFIdx extends RAF.SequenceFragment {

		public SequenceFragmentRAFIdx() {
			super();
		}

		public SequenceFragmentRAFIdx(int expectedLength) {
			super(expectedLength);
		}

		List<Integer> RAFindices = new ArrayList<Integer>();
		/**
		   appends legal (no . or ") character to sequence
		*/
		final public void append(char seqChar, int RAFidx) {
		    if ((seqChar != '.') &&
			(seqChar != '"')) {
			seqBuffer.append(seqChar);
			RAFindices.add(RAFidx);
			// System.out.println(RAFindices);
			if (seqChar == 'x')
			    xCount++;
		    }
		}

        /**
           appends another fragment
        */
        final public void append(SequenceFragmentRAFIdx f) {
            seqBuffer.append(f.seqBuffer);
            xCount += f.xCount;
	    RAFindices.addAll(f.RAFindices);
        }

	}


    /**
       returns sequence of part of a chain, according to source:
       <p/>
       1 = ATOM
       2 = SEQRES
       <p/>
       If residue ids are duplicated in the RAF (common with older
       PDB files), the range returned will be that between the FIRST
       instance of firstRes and the LAST instance of lastRes.
       Returns null if anything was not found or range is reversed.
    */
    final public static Pair<SequenceFragmentRAFIdx, SequenceFragmentRAFIdx> partialChainSeqPair(String body, String firstRes, String lastRes) {
        int firstResN = RAF.indexOf(body, firstRes, true);
        int lastResN = RAF.indexOf(body, lastRes, false);
        // System.out.println(firstRes+" "+firstResN);
        // System.out.println(lastRes+" "+lastResN);
        if ((firstResN == -1) || (lastResN == -1) || (firstResN > lastResN)) {
            return null;
        } else {
            return partialChainSeqPair(body, firstResN, lastResN);
        }      
    }

    /**
       returns sequence of part of a chain, according to source:
       <p/>
       1 = ATOM
       2 = SEQRES
       <p/>
       Residue ids are numeric, and 0-indexed.  lastRes must be
       greater than or equal to firstRes.  Indices are not checked;
       this will cause an exception if off either end.
    */
    final public static Pair<SequenceFragmentRAFIdx, SequenceFragmentRAFIdx> partialChainSeqPair(String body, 
                                                         int firstRes, 
                                                         int lastRes) throws IllegalArgumentException {
        SequenceFragmentRAFIdx rvATOM = new SequenceFragmentRAFIdx(lastRes - firstRes + 1);
        SequenceFragmentRAFIdx rvSEQRES = new SequenceFragmentRAFIdx(lastRes - firstRes + 1);
        for (int i = firstRes * 7; i <= lastRes * 7; i += 7) {
	    int RAFidx = Math.floorDiv(i, 7);
            rvATOM.append(body.charAt(i + 5), RAFidx);
            rvSEQRES.append(body.charAt(i + 6), RAFidx);
        }
        return new Pair<> (rvATOM, rvSEQRES);
    }

    /**
       returns sequence for a particular scop node, according
       to source:

       1 = ATOM
       2 = SEQRES

       and according to style:
       1 = single chain
       2 = multi-chain original style (OS)
       3 = multi-chain genetic domain (GD)

       If style 2, the order argument is a 0-based index saying which
       sequence to return.  For example, in "1avo A:,B:"
       order 0 would return 1avoA, and order 1 would return 1avoB.

       sequences are lower case, but separate parts of domains that
       are discontiguous on a chain (or on different chains, in the
       case of GD sequences) are separated by a capital X

       Returns null if error.
    */
    final public static Pair<SequenceFragmentRAFIdx, SequenceFragmentRAFIdx> domainSeq(int nodeID, int styleType, int order) {
        try {
            Statement stmt = LocalSQL.createStatement();
            ResultSet rs = stmt.executeQuery("select description from scop_node where id=" + nodeID);
            rs.next();
            String description = rs.getString(1).substring(5);
            rs.close();
            char chain = ' ';
            char lastChain = ' ';
            SequenceFragmentRAFIdx rvATOM = new SequenceFragmentRAFIdx();
            SequenceFragmentRAFIdx rvSEQRES = new SequenceFragmentRAFIdx();
            int currentChain = 0;
            boolean firstRegion = true;
            boolean addedFragment = false;
            String body = null;
            Pattern regionPattern = Pattern.compile("\\s*(\\S+?)-(\\S+)\\s*$");

            // System.out.println("node description is "+description);

            StringTokenizer st = new StringTokenizer(description, ",");
            while (st.hasMoreTokens()) {
                String region = st.nextToken();
                int pos = region.indexOf(':');
                if (pos == -1)
                    chain = ' ';
                else {
                    chain = region.charAt(pos - 1);
                    region = region.substring(pos + 1);
                }

                if (firstRegion) {
                    lastChain = chain;
                    firstRegion = false;
                }

                if (lastChain != chain)
                    currentChain++;

                if ((styleType != 2) || (currentChain == order)) {
                    // include this region
                    // System.out.println("using region "+region);
                    if ((body == null) || (chain != lastChain)) {
                        // get a new RAF body
                        System.out.println("getting RAF for "+nodeID);
			String query = "select raf_get_body(r.id) from raf r, pdb_chain c, link_pdb l, scop_node n where r.pdb_chain_id=c.id and l.pdb_chain_id=c.id and c.chain=\"" + chain + "\" and l.node_id=n.id and r.first_release_id is null and r.last_release_id is null and r.raf_version_id = 3 and n.id=" + nodeID + " group by r.line;";
                        rs = stmt.executeQuery(query);
                        if (!rs.next()) {
                            System.out.println("RAF not found");
                            rs = stmt.executeQuery("select raf_get_body(r.id) from raf r, pdb_chain c, link_pdb l, scop_node n where r.pdb_chain_id=c.id and l.pdb_chain_id=c.id and c.chain=\"" + Character.toLowerCase(chain) + "\" and l.node_id=n.id and n.release_id >= r.first_release_id and n.release_id <= r.last_release_id and n.id=" + nodeID);
                            if (rs.next()) {
                                System.out.println(nodeID + " " + description + " needs case fixed.");
                            } else {
                                rs.close();
                                stmt.close();
                                return null;
                            }
                        }
                        body = rs.getString(1);
			// System.out.println("body: "+ body);
                        if (rs.next()) {
                            System.out.println("error - raf logic wrong");
                            System.exit(1);
                        }
                        rs.close();
                    }
                    // figure out boundaries
                    String firstRes = null;
                    String lastRes = null;
                    Matcher m = regionPattern.matcher(region);
                    if (m.matches()) {
                        firstRes = m.group(1);
                        lastRes = m.group(2);
                    }
                    Pair<SequenceFragmentRAFIdx, SequenceFragmentRAFIdx> fPair = null;
		    SequenceFragmentRAFIdx fATOM = null;
		    SequenceFragmentRAFIdx fSEQRES = null;
		    

                    if ((firstRes != null) && (lastRes != null)) {
                        fPair = partialChainSeqPair(body, firstRes, lastRes);
			fATOM = fPair.getKey();
			fSEQRES = fPair.getValue();
		    } else {
                        fPair = wholeChainSeqPair(body);
			fATOM = fPair.getKey();
			fSEQRES = fPair.getValue();
                    }

                    // if error, blow the whole sequence
                    if (fPair == null) {
                        stmt.close();
                        return null;
                    }

                    // append f to end, with X if needed
                    if (addedFragment)
                        rvATOM.append('X');
                        rvSEQRES.append('X');
                    rvATOM.append(fATOM);
                    rvSEQRES.append(fSEQRES);
                    addedFragment = true;
                }

                lastChain = chain;
            }
            stmt.close();
            return new Pair<> (rvATOM, rvSEQRES);
        } catch (Exception e) {
            System.out.println("Exception: " + e.getMessage());
            e.printStackTrace();
        }
        return null;
    }

    /**
       returns sequence of entire chain, according to source:
       <p/>
       1 = ATOM
       2 = SEQRES
       3 = SEQRES, with first/last ATOM boundaries (pre-1.65)
       <p/>
       (as in astral_seq_source table)
    */
    final public static Pair<SequenceFragmentRAFIdx, SequenceFragmentRAFIdx> wholeChainSeqPair(String body) throws IllegalArgumentException {
        int l = body.length();
        SequenceFragmentRAFIdx rvATOM = new SequenceFragmentRAFIdx(l / 7);
        SequenceFragmentRAFIdx rvSEQRES = new SequenceFragmentRAFIdx(l / 7);
        for (int i = 0; i < l; i += 7) {
            char atomRes = body.charAt(i + 5);
            char seqRes = body.charAt(i + 6);
	    int RAFidx = Math.floorDiv(i, 7);

            if (atomRes != '.' && atomRes != '"')
                rvATOM.append(atomRes, RAFidx);

            if (seqRes != '.' && seqRes != '"')
                rvSEQRES.append(seqRes, RAFidx);
        }
        return new Pair<> (rvATOM, rvSEQRES);
    }


	public static void populateOneModel(int releaseID, int testModelId) throws Exception {


		    ArrayList<String> supfam = new ArrayList<String>();
		    Statement stmt = LocalSQL.createStatement();
		    String sqlquery = "select sccs from scop_node where release_id=" + releaseID + " and level_id=4 order by sccs";
		    ResultSet rs = stmt.executeQuery(sqlquery);
		    while (rs.next()) {
			    supfam.add(rs.getString(1));
		    }
		    rs.close();
		    ArrayList<String> allRep = new ArrayList<String>();
		    for (String s : supfam){
		           File su = new File("/h/slin/Rep"+releaseID+"/outf/"+s+"_rep.txt");
			   if (!su.exists()){
				   // System.out.println("no representatives found for " + s);
				   continue;
			   }
  			   Scanner myReader = new Scanner(su);
   			   while (myReader.hasNextLine()) {
				   String l = myReader.nextLine();
				   if (l.startsWith("[")) {
					   l=l.substring(1, l.length()-1).replace(" ", "");
					   String r = l.split(",")[0];
					   allRep.add(r);
			   }
			}
		}
		    String modelDir = getModelDir(testModelId);
		    for (String sid : allRep) {
			    populateOneComparison(modelDir, testModelId, sid, 15, releaseID);
		    }
        }

	private static String getModelDir(int modelId) throws Exception {
		    PreparedStatement stmt2 = LocalSQL.prepareStatement("select pdb_path from model_structure where id = ?;");
		    stmt2.setInt(1, modelId);
		    ResultSet modelRs = stmt2.executeQuery();
		    modelRs.next();
		    return modelRs.getString("pdb_path");
	}

	public static void testOneComparison(String sid) {
		int testModelId = 7615;
	        int releaseID = 17;
		try {
		    String modelDir = getModelDir(testModelId);
		    populateOneComparison(modelDir, testModelId, sid, 15, releaseID);
		} catch (Exception e) {
		    System.out.println("Exception: "+e.getMessage());
		    e.printStackTrace();
		}
	}
	public static void populateTest() {
		int testModelId = 7615;
	        int releaseID = 17;
		try {
		    populateOneModel(releaseID, testModelId);
		} catch (Exception e) {
		    System.out.println("Exception: "+e.getMessage());
		    e.printStackTrace();
		}
	}


	public static void main_gridengine(String[] args) {
		int testModelId = Integer.parseInt(args[0]);
	        int releaseID = 17;
		try {
		    populateOneModel(releaseID, testModelId);
		} catch (Exception e) {
		    System.out.println("Exception: "+e.getMessage());
		    e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		try {
		    LocalSQL.connectRW();
		    populateTest();
		} catch (Exception e) {
		    System.out.println("Exception: " + e.getMessage());
		    e.printStackTrace();
		}
	}

}
