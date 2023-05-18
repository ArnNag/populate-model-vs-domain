import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.util.Scanner;
import java.util.Arrays;
import java.util.*;

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

public class TestRAF {
	public static void main(String[] args) {
		int testDomainID = 830903;
		SequenceFragmentRAFIdx ATOM = domainSeq(testDomainID, 1, 1, 0);
		SequenceFragmentRAFIdx SEQRES = domainSeq(testDomainID, 2, 1, 0);
		System.out.println(ATOM.getSequence());
		System.out.println(ATOM.RAFindices);
		System.out.println(SEQRES.getSequence());
		System.out.println(SEQRES.RAFindices);
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
    final public static SequenceFragmentRAFIdx partialChainSeq(String body, int sourceType, String firstRes, String lastRes) {
        int firstResN = RAF.indexOf(body, firstRes, true);
        int lastResN = RAF.indexOf(body, lastRes, false);
        // System.out.println(firstRes+" "+firstResN);
        // System.out.println(lastRes+" "+lastResN);
        if ((firstResN == -1) || (lastResN == -1) || (firstResN > lastResN)) {
            return null;
        } else {
            return partialChainSeq(body, sourceType, firstResN, lastResN);
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
    final public static SequenceFragmentRAFIdx partialChainSeq(String body, 
                                                         int sourceType, 
                                                         int firstRes, 
                                                         int lastRes) throws IllegalArgumentException {
        SequenceFragmentRAFIdx rv = new SequenceFragmentRAFIdx(lastRes - firstRes + 1);
        for (int i = firstRes * 7; i <= lastRes * 7; i += 7) {
	    int RAFidx = Math.floorDiv(i, 7);
            if (sourceType == 1) {
                rv.append(body.charAt(i + 5), RAFidx);
            } else if (sourceType == 2) {
                rv.append(body.charAt(i + 6), RAFidx);
            } else {
                throw new IllegalArgumentException("'sourceType' must have value 1 or 2.");
            }
        }
        return rv;
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
    final public static SequenceFragmentRAFIdx domainSeq(int domainID, int sourceType, int styleType, int order) {
        try {
            Statement stmt = LocalSQL.createStatement();
            ResultSet rs = stmt.executeQuery("select description from scop_node where id=" + domainID);
            rs.next();
            String description = rs.getString(1).substring(5);
            rs.close();
            char chain = ' ';
            char lastChain = ' ';
            SequenceFragmentRAFIdx rv = new SequenceFragmentRAFIdx();
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
                        System.out.println("getting RAF for "+domainID);
			String query = "select raf_get_body(r.id) from raf r, pdb_chain c, link_pdb l, scop_node n where r.pdb_chain_id=c.id and l.pdb_chain_id=c.id and c.chain=\"" + chain + "\" and l.node_id=n.id and r.first_release_id is null and r.last_release_id is null and r.raf_version_id = 3 and n.id=" + domainID + ";";
			System.out.println(query);
                        rs = stmt.executeQuery(query);
                        if (!rs.next()) {
                            System.out.println("RAF not found");
                            rs = stmt.executeQuery("select raf_get_body(r.id) from raf r, pdb_chain c, link_pdb l, scop_node n where r.pdb_chain_id=c.id and l.pdb_chain_id=c.id and c.chain=\"" + Character.toLowerCase(chain) + "\" and l.node_id=n.id and n.release_id >= r.first_release_id and n.release_id <= r.last_release_id and n.id=" + domainID);
                            if (rs.next()) {
                                System.out.println(domainID + " " + description + " needs case fixed.");
                            } else {
                                rs.close();
                                stmt.close();
                                return null;
                            }
                        }
                        body = rs.getString(1);
			System.out.println("body: "+ body);
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
                    SequenceFragmentRAFIdx f = null;

                    if ((firstRes != null) && (lastRes != null)) {
                        f = partialChainSeq(body, sourceType, firstRes, lastRes);
		    } else {
                        // for SEQRES, we only want region within ATOMs
                        int st2 = sourceType;
                        if (st2 == 2) st2++;
                        f = wholeChainSeq(body, st2);
                    }

                    // if error, blow the whole sequence
                    if (f == null) {
                        stmt.close();
                        return null;
                    }

                    // append f to end, with X if needed
                    if (addedFragment)
                        rv.append('X');
                    rv.append(f);
                    addedFragment = true;
                }

                lastChain = chain;
            }
            stmt.close();
            return rv;
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
    final public static SequenceFragmentRAFIdx wholeChainSeq(String body, int sourceType) throws IllegalArgumentException {
        int l = body.length();
        SequenceFragmentRAFIdx rv = new SequenceFragmentRAFIdx(l / 7);
        boolean inSeq = true;
        int atomStart = 0;
        int atomEnd = l - 7;
        if (sourceType == 3) {
            for (int i = 0; i < l; i += 7) {
                if (body.charAt(i + 5) != '.') {
                    atomStart = i;
                    break;
                }
            }
            for (int i = l - 7; i >= 0; i -= 7) {
                if (body.charAt(i + 5) != '.') {
                    atomEnd = i;
                    break;
                }
            }
        }
        for (int i = 0; i < l; i += 7) {
            char atomRes = body.charAt(i + 5);
            char seqRes = body.charAt(i + 6);
	    int RAFidx = Math.floorDiv(i, 7);

            if (sourceType == 3) {
                if (i < atomStart || i > atomEnd) {
                    inSeq = false;
                } else {
                    inSeq = true;
                }
            }
            char seqChar;
            if (sourceType == 1) {
                seqChar = atomRes;
            } else if (sourceType == 2 || sourceType == 3) {
                seqChar = seqRes;
            } else {
                throw new IllegalArgumentException("'sourceType' must have value 1, 2, or 3.");
            }  
            if (inSeq && seqChar != '.' && seqChar != '"')
                rv.append(seqChar, RAFidx);
        }
        return rv;
    }

}
