from dash_dashboard.base_classes import BasePageLayout, BaseMain, BaseSideBar
import abc


class AliPageLayout(BasePageLayout, abc.ABC):
    top_bar_title = "LEED Viewer"


class AliMain(BaseMain, abc.ABC):
    pass


class AliSidebar(BaseSideBar, abc.ABC):
    pass