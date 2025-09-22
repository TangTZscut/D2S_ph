from pathlib import Path
import os
import shutil

class PathUtils:
    @staticmethod
    def get_current_path() -> str:
        """
        获取当前脚本文件所在的绝对路径
        :return: 字符串形式的路径
        """
        return str(Path(__file__).resolve().parent.parent)

    @staticmethod
    def make_dir_from_current(path: str) -> str:
        """
        在当前脚本文件所在的目录下创建目录
        :param path: 要创建的目录的相对路径
        :return: 字符串形式的绝对路径
        """
        current_path = PathUtils.get_current_path()
        abs_path = Path(current_path).joinpath(path).resolve()
        abs_path.mkdir(parents=True, exist_ok=True)
        return str(abs_path)
    
    @staticmethod
    def get_tmp_path(dir: str = "") -> str:
        tmp = os.path.join(PathUtils.get_current_path(), "tmp")
        if dir:
            tmp = os.path.join(tmp, dir)
        if not os.path.exists(tmp):
            os.makedirs(tmp)
        return tmp
    
    @staticmethod
    def del_tmp_dir() -> None:
        try:
            path = PathUtils.get_tmp_path()
            shutil.rmtree(path)
            print(f"已成功删除：{path}")
        except Exception as e:
            print(f"删除失败：{e}")