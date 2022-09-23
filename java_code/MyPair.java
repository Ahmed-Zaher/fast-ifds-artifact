/*
 * This class is mainly created to be used in other classes of the JTDec package. It is very similar to the C++ pair.
 * Please send bug reports and suggestions to goharshady@ist.ac.at
 */
public class MyPair<L,R> {
	private final L left;
	private final R right;
	
	public MyPair(L left, R right)
	{
		this.left=left;
		this.right=right;
	}
	
	public MyPair<R,L> Inverse()
	{
		MyPair<R,L> ret=new MyPair(right,left);
		return ret;
	}
	
	public L First()
	{
		return left;
	}
	public R Second()
	{
		return right;
	}
	
	public int hashCode()
	{
		return left.hashCode()^right.hashCode();
	}
	
	public boolean equals(Object obj)
	{
		if(!(obj instanceof MyPair))
			return false;
		MyPair given=(MyPair)obj;
		return this.left.equals(given.First()) && this.right.equals(given.Second());
	}
	
	public String toString()
	{
		return "("+left.toString()+", "+right.toString()+")";
	}

}