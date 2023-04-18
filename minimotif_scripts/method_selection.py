"""  """

def calc_mean_ic(ic):
    """ Calculate the mean IC score for the input matrix """
    shanon_sum = ic.sum(axis=1)
    # removes the first 4 numbers of the shannon IC.
    shanon_sum.drop(index=shanon_sum.index[:4], axis=0, inplace=True)
    # removes the last 3 numbers of the shannon ic
    shanon_sum.drop(shanon_sum.tail(3).index, inplace=True)
    return shanon_sum.mean()

