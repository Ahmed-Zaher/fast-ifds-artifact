import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import soot.*;
import soot.jimple.*;
import soot.jimple.toolkits.callgraph.CallGraph;
import soot.jimple.toolkits.callgraph.Edge;
import soot.options.Options;
import soot.toolkits.graph.BriefUnitGraph;

/*
	Input:

	A DaCapo benchmark

	Purpose:

	Get the necessary details that make up the corresponding IFDS Arena


	Output: a supergraph of the following format:
	{
		n_G m_G n_H // n_G is number of CFG nodes which well label from 0 to n_G-1; m_G is the number of same-CFG edges;
					// and n_H is number of procedures which we'll label from 0 to n_G-1
		procOf[0] procOf[1] .. procOf[n_G-1] // n_G space-separated numbers, where procOf[i] is in [0, n_H) which denotes
											 // the label of procedure that CFG node i lies in
		u_1 v_1
		u_2 v_2
		   .
		   .
		u_{m_G} v_{m_G} // m_G edges specifying the same-CFG edges, has the form "u v" for u and v in [0, n_G)
		m_H		// number of call-graph edges
		u_1 p_1
		u_2 p_2
		   .
		   .
		u_{m_H} p_{m_H} // m_H edges specifying the procedure-call edges, has the form "u p"
						// for u in [0, n_G) and p in [0, n_H), meaning that CFG node u calls procedure p
		D[0] D[1] .. D[n_H-1] // D[p] denoting the domain size for procedure p
		// the following is omitted in case of reachability analysis
		p |params[0]| params[0]
		p |params[1]| params[1]
		    .
		    .
		p |params[n_H-1]| params[n_H-1] // params[p] are the indices of local variables that are in the domain of p
										// (their values are in [1, D[p]]), and are passed as parameters by a callee procedure
		c |calls[0]| calls[0]
		c |calls[1]| calls[1]
		    .
		    .
		c |calls[n_G-1]| calls[n_G-1] // where calls[u] is a list of pairs (x, y) meaning that we pass local variable x (w.r.t. procOf[u])
									  // to be assigned to local variable y (w.r.t called procedure)
		rets[0] rets[1] .. rets[n_G-1] // if u is a return statement then ret[u] = id of local variable returned if the return value 
									   // is indeed a local variable and in the domain; otherwise ret[u] = -1
		a |assignments[0]| assignments[0]
		a |assignments[1]| assignments[1]
		    .
		    .
		a |assignments[n_G-1]| assignments[n_G-1] // if u is an assignment statement and has the form x := .. where x is a local variable
												  // in the domain, then assignment[u] is [x, y1, y2, ..] where y's are local variables in 
												  // domain that appear in RHS; otherwise assignments[u] is empty
	}

 */


public class FastIFDS extends SceneTransformer {

	static String classPath, benchmarkName;

	public static void main(String[] args) {

		benchmarkName = args[0];

		System.out.println(args[0]);

		if (args[0].equals("avrora")) {
			args[0] = "avrora.Main";
			classPath = System.getProperty("user.dir") + File.separator + "lib" + File.separator + "dacapo-9.12-MR1-bach" + File.separator + "jar" + File.separator + "avrora-cvs-20091224.jar";
		} else if (args[0].equals("batik")) {
			args[0] = "org.apache.batik.apps.svgbrowser.Main";
			classPath = System.getProperty("user.dir") + File.separator + "lib" + File.separator + "dacapo-9.12-MR1-bach" + File.separator + "jar" + File.separator + "batik-all.jar";
		} else if (args[0].equals("h2")) {
			args[0] = "org.h2.tools.Console";
			classPath = System.getProperty("user.dir") + File.separator + "lib" + File.separator + "dacapo-9.12-MR1-bach" + File.separator + "jar" + File.separator + "h2-1.2.121.jar";
		} else if (args[0].equals("sunflow")) {
			args[0] = "SunflowGUI";
			classPath = System.getProperty("user.dir") + File.separator + "lib" + File.separator + "dacapo-9.12-MR1-bach" + File.separator + "jar" + File.separator + "sunflow-0.07.2.jar";
		} else {
			classPath = System.getProperty("user.dir") + File.separator + "lib" + File.separator + "dacapo-2006-10-MR2.jar";
			args[0] = "dacapo." + args[0] + ".Main";
		}

		System.out.println(args[0]);


		Options.v().set_soot_classpath(classPath);


		// Enable whole-program mode
		Options.v().set_whole_program(true);
		Options.v().set_app(true);
		Options.v().set_prepend_classpath(true);
		Options.v().set_allow_phantom_refs(true);


		// Call-graph options
		Options.v().setPhaseOption("cg.cha","enabled:false");

		// Enable SPARK call-graph construction
		Options.v().setPhaseOption("cg.spark","enabled:true");
		Options.v().setPhaseOption("cg.spark","verbose:true");
		Options.v().setPhaseOption("cg.spark","on-fly-cg:true");



		// Load the main class
		SootClass c = Scene.v().loadClass(args[0], SootClass.BODIES);
		c.setApplicationClass();



		// Load the "main" method of the main class and set it as a Soot entry point

		SootMethod entryPoint = c.getMethodByName("main");
		List<SootMethod> entryPoints = new ArrayList<SootMethod>();
		entryPoints.add(entryPoint);
		Scene.v().setEntryPoints(entryPoints);


		PackManager.v().getPack("wjtp").add(new Transform("wjtp.main", new FastIFDS()));

		soot.Main.main(args);

	}

	@Override
	protected void internalTransform(String phaseName, Map<String, String> options) {

		Scene.v().loadBasicClasses();


		try {
			run();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}



	public void run() throws IOException, InterruptedException {


		List<SootMethod> applicationMethods = soot.EntryPoints.v().methodsOfApplicationClasses();


		// Getting call-graph
		CallGraph cg = Scene.v().getCallGraph();


		HashSet<MyPair<Integer, Integer>> edgeListH = new HashSet<MyPair<Integer, Integer>>();

		ArrayList<MyPair<Integer, Integer>> edgeListG = new ArrayList<MyPair<Integer, Integer>> ();

		Integer n_H = 0, n_G = 0;
		HashMap<SootMethod, Integer> mappingProc = new HashMap<SootMethod, Integer>();
		HashMap<Unit, Integer> mappingUnit = new HashMap<Unit, Integer>();
		HashMap<Integer, Integer> procOf = new HashMap<Integer, Integer>();



		// Finding n_H, n_G, procOf[.]
		for (SootMethod method : applicationMethods) {
			mappingProc.put(method, n_H);
			n_H = n_H + 1;

			BriefUnitGraph cfg = new BriefUnitGraph(method.retrieveActiveBody());
			Iterator<Unit> iter = cfg.iterator();

			while (iter.hasNext()) {
				Unit u = (Unit) iter.next();
				mappingUnit.put(u, n_G);
				procOf.put(mappingUnit.get(u), mappingProc.get(method));
				n_G = n_G + 1;
			}
		}


		// Building G and H


		// G


		for (SootMethod method : applicationMethods) {

			BriefUnitGraph cfg = new BriefUnitGraph(method.retrieveActiveBody());
			Iterator<Unit> iter = cfg.iterator();

			while (iter.hasNext()) {
				Unit u = (Unit) iter.next();
				List<Unit> succ = cfg.getSuccsOf(u);

				for (int j = 0; j < succ.size(); ++j) {
					Unit v = succ.get(j);
					assert mappingUnit.containsKey(u);
					assert mappingUnit.containsKey(v);
					edgeListG.add(new MyPair<Integer, Integer>(mappingUnit.get(u), mappingUnit.get(v)));
				}

			}
		}

		// H

		for (SootMethod method : applicationMethods) {

			for (Iterator<Edge> it = cg.edgesOutOf(method); it.hasNext(); ) {

				Edge edge = it.next();

				if (!mappingProc.containsKey(edge.tgt()))
					continue;

				if (edge.srcStmt() != null && edge.srcStmt().containsInvokeExpr() && edge.tgt().equals(edge.srcStmt().getInvokeExpr().getMethod())) {
					assert mappingUnit.containsKey(edge.srcStmt());
					edgeListH.add(new MyPair<Integer, Integer>(mappingUnit.get(edge.srcStmt()), mappingProc.get(edge.tgt())));
				}
			}

		}

		System.out.println("n_G = " + n_G + ", m_G = " + edgeListG.size());
		System.out.println("n_H = " + n_H + ", m_H = " + edgeListH.size());


		// Doing the analysis, finding D and flows

		List<String> analyses = new ArrayList<String>();

		analyses.add("uninit_var");
		analyses.add("reachability");
		analyses.add("null_ptr");

		// inDeg[mappingProc[p]] = in-degree of procedure p in H, same with outDeg
		HashMap<Integer, Integer> inDegProc = new HashMap<>();
		HashMap<Integer, Integer> outDegProc = new HashMap<>();

		for (SootMethod method : applicationMethods) {
			inDegProc.put(mappingProc.get(method), 0);
			outDegProc.put(mappingProc.get(method), 0);
		}

		for (MyPair<Integer, Integer> pr : edgeListH) {
			inDegProc.put(pr.Second(), inDegProc.get(pr.Second()) + 1);
			outDegProc.put(procOf.get(pr.First()), outDegProc.get(procOf.get(pr.First())) + 1);
		}


		for (String analysis : analyses) {


			Integer badProcCnt = 0;

			PrintWriter writer = new PrintWriter("../java_output/" + analysis + "/" + benchmarkName + ".txt", "UTF-8");

			System.out.println("Doing analysis: " + analysis);

			HashMap<Integer, Integer> D = new HashMap<Integer, Integer>();

			// dom[p] is a mapping from v to i meaning that value v maps to i in [1, dom[p]]
			HashMap<Integer, HashMap<Value, Integer>> dom = new HashMap<Integer, HashMap<Value, Integer>>();

			// params[p] is a mapping from i to j meaning that i'th argument maps to a local variable in dom[p] and has index j w.r.t. dom[p]
			// the index of an argument is not apparent in output, only used internally
			HashMap<Integer, HashSet<Integer>> params = new HashMap<Integer, HashSet<Integer>>();
			HashMap<Integer, HashMap<Integer, Integer>> calls = new HashMap<Integer, HashMap<Integer, Integer>>();
			HashMap<Integer, Integer> rets = new HashMap<Integer, Integer>();
			HashMap<Integer, ArrayList<Integer>> assignments = new HashMap<Integer, ArrayList<Integer>>();


			// reachability analysis
			if (analysis.equals("reachability")) {
				for (int i = 0; i < n_H; i++)
					D.put(i, 1);
			} else {

				// building dom and D
				for (SootMethod method : applicationMethods) {
					Body b = method.retrieveActiveBody();
					BriefUnitGraph cfg = new BriefUnitGraph(b);
					Iterator<Unit> iter = cfg.iterator();
					Integer p = mappingProc.get(method);
					dom.put(p, new HashMap<Value, Integer>());
					Integer domSz = 0;

					while (iter.hasNext()) {
						Unit u = (Unit) iter.next();
						List<ValueBox> L = u.getUseAndDefBoxes();
						for (ValueBox vb : L) {
							Value val = vb.getValue();
							if (val instanceof Local) {

								boolean inDom = false;

								if (analysis.equals("uninit_var"))
									inDom = true;
								if (analysis.equals("null_ptr"))
									inDom = (val.getType() instanceof RefLikeType);

								if (inDom && !dom.get(p).containsKey(val)) {
									domSz++;
									dom.get(p).put(val, domSz);
								}
							}
						}
					}

					D.put(p, domSz);

				}

				// params
				for (SootMethod method : applicationMethods) {
					Body b = method.retrieveActiveBody();
					BriefUnitGraph cfg = new BriefUnitGraph(b);
					Integer p = mappingProc.get(method);
					params.put(p, new HashSet<Integer>());


					for (int i = 0; i < method.getParameterCount(); ++i) {
						Local lcl = method.getActiveBody().getParameterLocal(i);
						if (dom.get(p).containsKey(lcl)) {
							params.get(p).add(dom.get(p).get(lcl));
						}
					}
				}


				// for each call, record which parameters are passed
				for (SootMethod method : applicationMethods) {
					Body b = method.retrieveActiveBody();
					BriefUnitGraph cfg = new BriefUnitGraph(b);
					Iterator<Unit> iter = cfg.iterator();
					Integer p = mappingProc.get(method);

					while (iter.hasNext()) {
						Unit u = (Unit) iter.next();
						Stmt stmt = (Stmt) u;
						Integer uMap = mappingUnit.get(u);

						calls.put(uMap, new HashMap<Integer, Integer>());


						if (stmt.containsInvokeExpr() && mappingProc.containsKey(stmt.getInvokeExpr().getMethod())) {
							InvokeExpr invokeExpr = stmt.getInvokeExpr();
							SootMethod destinationMethod = invokeExpr.getMethod();
							Integer pp = mappingProc.get(destinationMethod);
							assert mappingUnit.containsKey(stmt);

							if (edgeListH.contains(new MyPair<Integer, Integer>(mappingUnit.get(stmt), pp))) {
								List<Value> args = invokeExpr.getArgs();
								List<Local> localArguments = new ArrayList();
								Iterator iter2 = args.iterator();

								while (iter2.hasNext()) {
									Value value = (Value) iter2.next();
									if (value instanceof Local) {
										localArguments.add((Local)value);
									}
								}

								for (int i = 0; i < localArguments.size(); i++) {
									Local lcl = localArguments.get(i);
									Local correspondingLcl = destinationMethod.getActiveBody().getParameterLocal(i);
									if (dom.get(p).containsKey(lcl) && dom.get(pp).containsKey(correspondingLcl)) {
										calls.get(uMap).put(dom.get(p).get(lcl), dom.get(pp).get(correspondingLcl));
									}
								}
							}
						}
					}

				}

				// for each return, record which parameter is returned w.r.t. callee procedure
				for (SootMethod method : applicationMethods) {
					Body b = method.retrieveActiveBody();
					BriefUnitGraph cfg = new BriefUnitGraph(b);
					Iterator<Unit> iter = cfg.iterator();
					Integer p = mappingProc.get(method);


					while (iter.hasNext()) {
						Unit u = (Unit) iter.next();
						Integer uMap = mappingUnit.get(u);

						rets.put(uMap, -1);

						if (u instanceof ReturnStmt) {
							ReturnStmt retStmt = (ReturnStmt) u;
							if (retStmt.getOp() instanceof Local) {
								if (dom.get(p).containsKey(retStmt.getOp())) {
									rets.put(uMap, dom.get(p).get(retStmt.getOp()));
								}
							}
						}
					}
				}

				// for each assignment, record the label of LHS as well as labels of local variables appearing in RHS
				for (SootMethod method : applicationMethods) {
					Body b = method.retrieveActiveBody();
					BriefUnitGraph cfg = new BriefUnitGraph(b);
					Iterator<Unit> iter = cfg.iterator();
					Integer p = mappingProc.get(method);


					while (iter.hasNext()) {
						Unit u = (Unit) iter.next();
						Integer uMap = mappingUnit.get(u);

						assignments.put(uMap, new ArrayList<Integer>());

						if (u instanceof DefinitionStmt) {
							DefinitionStmt definition = (DefinitionStmt) u;
							Value leftOp = definition.getLeftOp();
							if (leftOp instanceof Local) {
								Local leftOpLocal = (Local) leftOp;
								if (dom.get(p).containsKey(leftOpLocal)) {
									assignments.get(uMap).add(dom.get(p).get(leftOpLocal));
									List<ValueBox> useBoxes = definition.getUseBoxes();
									for (ValueBox vb : useBoxes) {
										if (dom.get(p).containsKey(vb.getValue())) {
											assignments.get(uMap).add(dom.get(p).get(vb.getValue()));
										}
									}
								}
							}
						}

					}
				}
			}



			// Printing G, H

			writer.println(n_G + " " + edgeListG.size() + " " + n_H);

			for (int u = 0; u < n_G; u++) {
				writer.print(procOf.get(u));
				if (u + 1 < n_G)
					writer.print(" ");
			}
			writer.println();

			for (MyPair<Integer, Integer> pr : edgeListG)
				writer.println(pr.First() + " " + pr.Second());

			writer.println(edgeListH.size());

			for (MyPair<Integer, Integer> pr : edgeListH)
				writer.println(pr.First() + " " + pr.Second());



			// Printing D and flow information


			writer.print(D.get(0));
			for (int i = 1; i < n_H; i++) {
				writer.print(" " + D.get(i));
			}
			writer.println();

			if (!analysis.equals("reachability")) {
				for (int p = 0; p < n_H; p++) {
					writer.print("p ");
					writer.print(params.get(p).size());
					for (Integer param : params.get(p)) {
						writer.print(" " + param);
					}
					writer.println();
				}

				for (int i = 0; i < n_G; i++) {
					writer.print("c ");
					writer.print(calls.get(i).size());
					for (Map.Entry<Integer, Integer> entry : calls.get(i).entrySet()) {
						Integer x = entry.getKey();
						Integer y = entry.getValue();
						writer.print(" " + x + " " + y);
					}
					writer.println();
				}


				writer.print(rets.get(0));
				for (int i = 1; i < n_G; i++) {
					assert rets.containsKey(i);
					writer.print(" " + rets.get(i));
				}
				writer.println();


				for (int i = 0; i < n_G; i++) {
					writer.print("a ");
					writer.print(assignments.get(i).size());
					for (Integer value : assignments.get(i)) {
						writer.print(" " + value);
					}
					writer.println();
				}
			}
			writer.println("done");
			writer.close();
		}
	}
}
